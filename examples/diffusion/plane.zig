const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const cg = flux.math.cg;
const operator_context_mod = flux.operators.context;
const evolution_mod = flux.integrators.evolution;

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
    reduced_matrix: sparse.CsrMatrix(f64),
    reduced_index: []u32,
    interior_vertices: []u32,
    boundary_mask: []bool,
    diagonal: []f64,
    scratch: cg.Scratch,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const Mesh2D,
        operator_context: *operator_context_mod.OperatorContext(Mesh2D),
        dt: f64,
    ) !HeatSystem {
        const laplacian = try operator_context.laplacian(0);
        const stiffness = laplacian.stiffness;
        const masses = mesh.vertices.slice().items(.dual_volume);

        const boundary_mask = try boundaryVertexMask(allocator, mesh);
        errdefer allocator.free(boundary_mask);

        const reduced_index = try allocator.alloc(u32, mesh.num_vertices());
        errdefer allocator.free(reduced_index);
        @memset(reduced_index, std.math.maxInt(u32));

        var interior_count: u32 = 0;
        for (0..mesh.num_vertices()) |vertex_idx_usize| {
            const vertex_idx: u32 = @intCast(vertex_idx_usize);
            if (boundary_mask[vertex_idx]) continue;
            reduced_index[vertex_idx] = interior_count;
            interior_count += 1;
        }

        const interior_vertices = try allocator.alloc(u32, interior_count);
        errdefer allocator.free(interior_vertices);
        for (0..mesh.num_vertices()) |vertex_idx_usize| {
            const vertex_idx: u32 = @intCast(vertex_idx_usize);
            if (boundary_mask[vertex_idx]) continue;
            interior_vertices[reduced_index[vertex_idx]] = vertex_idx;
        }

        var triplets = sparse.TripletAssembler(f64).init(interior_count, interior_count);
        defer triplets.deinit(allocator);
        for (0..mesh.num_vertices()) |row_idx_usize| {
            const row_idx: u32 = @intCast(row_idx_usize);
            if (boundary_mask[row_idx]) continue;

            const reduced_row = reduced_index[row_idx];
            try triplets.addEntry(allocator, reduced_row, reduced_row, masses[row_idx]);

            const row = stiffness.row(row_idx);
            for (row.cols, row.vals) |col_idx, value| {
                if (boundary_mask[col_idx]) continue;
                try triplets.addEntry(allocator, reduced_row, reduced_index[col_idx], dt * value);
            }
        }

        var reduced_matrix = try triplets.build(allocator);
        errdefer reduced_matrix.deinit(allocator);

        const diagonal = try allocator.alloc(f64, interior_count);
        errdefer allocator.free(diagonal);
        for (0..interior_count) |row_idx_usize| {
            const row_idx: u32 = @intCast(row_idx_usize);
            diagonal[row_idx] = diagonalEntry(reduced_matrix, row_idx);
            std.debug.assert(diagonal[row_idx] > 0.0);
        }

        var scratch = try cg.Scratch.init(allocator, interior_count);
        errdefer scratch.deinit(allocator);

        return .{
            .reduced_matrix = reduced_matrix,
            .reduced_index = reduced_index,
            .interior_vertices = interior_vertices,
            .boundary_mask = boundary_mask,
            .diagonal = diagonal,
            .scratch = scratch,
        };
    }

    pub fn deinit(self: *HeatSystem, allocator: std.mem.Allocator) void {
        self.scratch.deinit(allocator);
        allocator.free(self.diagonal);
        allocator.free(self.boundary_mask);
        allocator.free(self.interior_vertices);
        allocator.free(self.reduced_index);
        self.reduced_matrix.deinit(allocator);
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

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    const operator_context = try operator_context_mod.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();

    const total_time = final_time_override orelse (@as(f64, @floatFromInt(config.steps)) * config.dt());
    const dt = total_time / @as(f64, @floatFromInt(config.steps));
    var heat_system = try HeatSystem.init(allocator, &mesh, operator_context, dt);
    defer heat_system.deinit(allocator);

    var state = try VertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeState(&mesh, state.values, initial_condition, 0.0);

    const Evolution = evolution_mod.ExactEvolution(
        Mesh2D,
        f64,
        HeatStepperBuilder,
        HeatExactInitializer,
        HeatErrorMeasure,
    );
    var evolution = try Evolution.init(
        allocator,
        &mesh,
        state.values,
        HeatStepperBuilder{
            .mesh = &mesh,
            .heat_system = &heat_system,
        },
        HeatExactInitializer{ .initial_condition = initial_condition },
        HeatErrorMeasure{},
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

fn boundaryVertexMask(allocator: std.mem.Allocator, mesh: *const Mesh2D) ![]bool {
    const mask = try allocator.alloc(bool, mesh.num_vertices());
    @memset(mask, false);

    const edge_vertices = mesh.simplices(1).items(.vertices);
    for (mesh.boundary_edges) |edge_idx| {
        const edge = edge_vertices[edge_idx];
        mask[edge[0]] = true;
        mask[edge[1]] = true;
    }
    return mask;
}

fn diagonalEntry(matrix: sparse.CsrMatrix(f64), row_idx: u32) f64 {
    const row = matrix.row(row_idx);
    for (row.cols, row.vals) |col_idx, value| {
        if (col_idx == row_idx) return value;
    }
    unreachable;
}

fn seedReducedSolution(
    heat_system: *const HeatSystem,
    full_state: []const f64,
    reduced_solution: []f64,
) void {
    std.debug.assert(reduced_solution.len == heat_system.interior_vertices.len);
    for (heat_system.interior_vertices, 0..) |vertex_idx, interior_idx| {
        reduced_solution[interior_idx] = full_state[vertex_idx];
    }
}

const HeatStepperBuilder = struct {
    pub const Stepper = HeatStepper;

    mesh: *const Mesh2D,
    heat_system: *const HeatSystem,

    pub fn initStepper(self: @This(), allocator: std.mem.Allocator, state_values: []f64) !Stepper {
        const len = self.heat_system.interior_vertices.len;
        const stepper = Stepper{
            .mesh = self.mesh,
            .heat_system = self.heat_system,
            .state_values = state_values,
            .reduced_rhs = try allocator.alloc(f64, len),
            .reduced_solution = try allocator.alloc(f64, len),
        };
        seedReducedSolution(self.heat_system, state_values, stepper.reduced_solution);
        return stepper;
    }
};

const HeatStepper = struct {
    mesh: *const Mesh2D,
    heat_system: *const HeatSystem,
    state_values: []f64,
    reduced_rhs: []f64,
    reduced_solution: []f64,

    pub fn step(self: *@This(), allocator: std.mem.Allocator) !void {
        try stepBackwardEuler(
            allocator,
            self.mesh,
            self.heat_system,
            self.state_values,
            self.reduced_rhs,
            self.reduced_solution,
        );
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.reduced_solution);
        allocator.free(self.reduced_rhs);
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

fn stepBackwardEuler(
    allocator: std.mem.Allocator,
    mesh: *const Mesh2D,
    heat_system: *const HeatSystem,
    state_values: []f64,
    reduced_rhs: []f64,
    reduced_solution: []f64,
) !void {
    _ = allocator;
    const masses = mesh.vertices.slice().items(.dual_volume);
    std.debug.assert(reduced_rhs.len == heat_system.interior_vertices.len);
    std.debug.assert(reduced_solution.len == heat_system.interior_vertices.len);

    for (heat_system.interior_vertices, 0..) |vertex_idx, interior_idx| {
        reduced_rhs[interior_idx] = masses[vertex_idx] * state_values[vertex_idx];
    }

    var preconditioner = cg.DiagonalPreconditioner{ .diagonal = heat_system.diagonal };
    const solve_result = cg.solve(
        heat_system.reduced_matrix,
        reduced_rhs,
        reduced_solution,
        1e-10,
        4000,
        &preconditioner,
        heat_system.scratch,
    );
    if (!solve_result.converged) return error.ConjugateGradientDidNotConverge;

    for (heat_system.boundary_mask, 0..) |is_boundary, vertex_idx| {
        if (is_boundary) state_values[vertex_idx] = 0.0;
    }
    for (heat_system.interior_vertices, 0..) |vertex_idx, interior_idx| {
        state_values[vertex_idx] = reduced_solution[interior_idx];
    }
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
