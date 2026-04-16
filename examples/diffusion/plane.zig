const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const linear_system = flux.math.linear_system;
const feec_context_mod = flux.operators.feec.context;
const evolution_mod = flux.evolution;

pub const Mesh2D = flux.topology.Mesh(2, 2);
pub const VertexField = flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal);

const convergence_time = 0.02;

const InitialCondition = enum {
    zero,
    sine_mode,
};

pub const ConfigImpl = struct {
    grid: u32 = 8,
    steps: u32 = 8,
    domain: f64 = 1.0,
    dt_scale: f64 = 0.1,
    dt_override: ?f64 = null,
    output_dir: []const u8 = "output/heat",
    frames: u32 = 4,

    pub fn dt(self: ConfigImpl) f64 {
        if (self.dt_override) |value| return value;
        const h = self.domain / @as(f64, @floatFromInt(self.grid));
        return self.dt_scale * h * h;
    }
};

pub const RunResultImpl = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

pub const ConvergenceResultImpl = struct {
    grid: u32,
    l2_error: f64,
};

const HeatSystem = struct {
    linear_system: linear_system.LinearSystem,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const Mesh2D,
        operator_context: *feec_context_mod.OperatorContext(Mesh2D),
        dt: f64,
    ) !HeatSystem {
        const laplacian = try operator_context.laplacian(0);
        const stiffness = laplacian.stiffness;
        const masses = mesh.vertices.slice().items(.dual_volume);

        const elimination_map = try linear_system.EliminationMap.initBoundary(Mesh2D, allocator, mesh, 0);

        var triplets = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
        defer triplets.deinit(allocator);
        for (0..mesh.num_vertices()) |row_idx_usize| {
            const row_idx: u32 = @intCast(row_idx_usize);
            try triplets.addEntry(allocator, row_idx, row_idx, masses[row_idx]);

            const row = stiffness.row(row_idx);
            for (row.cols, row.vals) |col_idx, value| {
                try triplets.addEntry(allocator, row_idx, col_idx, dt * value);
            }
        }

        var full_matrix = try triplets.build(allocator);
        errdefer full_matrix.deinit(allocator);

        var system_runtime = try linear_system.LinearSystem.eliminate(
            allocator,
            full_matrix,
            elimination_map,
            .{},
        );
        errdefer system_runtime.deinit(allocator);
        full_matrix.deinit(allocator);

        return .{
            .linear_system = system_runtime,
        };
    }

    pub fn deinit(self: *HeatSystem, allocator: std.mem.Allocator) void {
        self.linear_system.deinit(allocator);
    }
};

pub fn runImpl(
    allocator: std.mem.Allocator,
    config: ConfigImpl,
    writer: anytype,
) !RunResultImpl {
    const result = try simulateCase(allocator, config, .sine_mode, null, writer);
    try writer.print(
        "heat: grid={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.grid, config.steps, config.dt(), result.l2_error },
    );
    return result;
}

pub fn runConvergenceStudyImpl(
    allocator: std.mem.Allocator,
    grids: []const u32,
) ![]ConvergenceResultImpl {
    const Runner = struct {
        pub fn run(allocator_inner: std.mem.Allocator, grid: u32) !ConvergenceResultImpl {
            const config = convergenceConfig(grid);
            const run_result = try simulateCase(allocator_inner, config, .sine_mode, convergence_time, null);
            return .{
                .grid = grid,
                .l2_error = run_result.l2_error,
            };
        }
    };

    return common.runConvergenceStudy(ConvergenceResultImpl, u32, allocator, grids, Runner{});
}

fn simulateCase(
    allocator: std.mem.Allocator,
    config: ConfigImpl,
    initial_condition: InitialCondition,
    final_time_override: ?f64,
    progress_writer: ?*std.Io.Writer,
) !RunResultImpl {
    std.debug.assert(config.grid > 0);
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.domain > 0.0);
    std.debug.assert(config.dt_scale > 0.0);

    var mesh = try Mesh2D.plane(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    const operator_context = try feec_context_mod.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();

    const total_time = final_time_override orelse (@as(f64, @floatFromInt(config.steps)) * config.dt());
    const dt = total_time / @as(f64, @floatFromInt(config.steps));
    var heat_system = try HeatSystem.init(allocator, &mesh, operator_context, dt);
    defer heat_system.deinit(allocator);

    var state = try VertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeState(&mesh, state.values, initial_condition, 0.0);

    const stepper = HeatStepper.init(&mesh, &heat_system, state.values);
    const aux = try HeatReferenceAux.init(
        allocator,
        &mesh,
        state.values.len,
        HeatExactInitializer{ .initial_condition = initial_condition },
        HeatErrorMeasure{},
    );

    const Evolution = evolution_mod.Evolution([]f64, HeatStepper, HeatReferenceAux);
    var evolution = Evolution.init(
        allocator,
        state.values,
        stepper,
        aux,
    );
    defer evolution.deinit();

    const loop_result = try common.runExactEvolutionLoop(
        Mesh2D,
        allocator,
        &mesh,
        &evolution,
        .{
            .steps = config.steps,
            .dt = dt,
            .final_time = total_time,
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = "heat",
            .progress_writer = progress_writer,
        },
    );

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .steps = config.steps,
        .snapshot_count = loop_result.snapshot_count,
        .l2_error = evolution.l2Error(),
    };
}

fn convergenceConfig(grid: u32) ConfigImpl {
    const probe = ConfigImpl{ .grid = grid };
    const dt_max = probe.dt();
    const step_count = @max(1, @as(u32, @intFromFloat(@ceil(convergence_time / dt_max))));
    return .{
        .grid = grid,
        .steps = step_count,
        .domain = 1.0,
        .dt_scale = probe.dt_scale,
        .frames = 0,
        .output_dir = "output/heat-convergence",
    };
}

const HeatStepper = struct {
    mesh: *const Mesh2D,
    heat_system: *HeatSystem,
    state_values: []f64,

    pub fn init(mesh: *const Mesh2D, heat_system: *HeatSystem, state_values: []f64) @This() {
        heat_system.linear_system.seedSolutionFromFull(state_values);
        return .{
            .mesh = mesh,
            .heat_system = heat_system,
            .state_values = state_values,
        };
    }

    pub fn step(self: *@This(), allocator: std.mem.Allocator) !void {
        try stepBackwardEuler(
            allocator,
            self.mesh,
            self.heat_system,
            self.state_values,
        );
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        _ = self;
        _ = allocator;
    }
};

const HeatErrorMeasure = struct {
    pub fn compute(_: @This(), mesh: *const Mesh2D, approx: []const f64, exact: []const f64) f64 {
        return weightedL2Error(mesh, approx, exact);
    }
};

const HeatExactInitializer = struct {
    initial_condition: InitialCondition,

    pub fn fill(self: @This(), mesh: *const Mesh2D, values: []f64, time: f64) void {
        initializeState(mesh, values, self.initial_condition, time);
    }
};
const HeatReferenceAux = evolution_mod.ReferenceAux(Mesh2D, HeatExactInitializer, HeatErrorMeasure);

fn stepBackwardEuler(
    allocator: std.mem.Allocator,
    mesh: *const Mesh2D,
    heat_system: *HeatSystem,
    state_values: []f64,
) !void {
    _ = allocator;
    const masses = mesh.vertices.slice().items(.dual_volume);
    const full_rhs = heat_system.linear_system.fullRhsValues();
    std.debug.assert(full_rhs.len == state_values.len);

    for (full_rhs, masses, state_values) |*rhs_value, mass, state_value| {
        rhs_value.* = mass * state_value;
    }

    _ = try heat_system.linear_system.solveHomogeneous(state_values);
}

fn initializeState(
    mesh: *const Mesh2D,
    values: []f64,
    initial_condition: InitialCondition,
    time: f64,
) void {
    const coords = mesh.vertices.slice().items(.coords);
    std.debug.assert(values.len == coords.len);
    for (coords, values) |coord, *value| {
        value.* = switch (initial_condition) {
            .zero => 0.0,
            .sine_mode => exactSolution(coord, time),
        };
    }
}

fn exactSolution(coords: [2]f64, time: f64) f64 {
    const eigenvalue = 2.0 * std.math.pi * std.math.pi;
    return std.math.exp(-eigenvalue * time) *
        std.math.sin(std.math.pi * coords[0]) *
        std.math.sin(std.math.pi * coords[1]);
}

fn weightedL2Error(mesh: *const Mesh2D, approx: []const f64, exact: []const f64) f64 {
    std.debug.assert(approx.len == exact.len);
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    var error_sq: f64 = 0.0;
    var measure: f64 = 0.0;
    for (approx, exact, dual_volumes) |approx_value, exact_value, dual_volume| {
        const diff = approx_value - exact_value;
        error_sq += diff * diff * dual_volume;
        measure += dual_volume;
    }
    return @sqrt(error_sq / measure);
}
