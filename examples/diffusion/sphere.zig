const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const linear_system = flux.math.linear_system;
const feec_context_mod = flux.operators.feec.context;
const evolution_mod = flux.evolution;

pub const SurfaceMesh = flux.topology.Mesh(3, 2);
pub const VertexField = flux.forms.Cochain(SurfaceMesh, 0, flux.forms.Primal);
const SurfaceOperatorContext = feec_context_mod.OperatorContext(SurfaceMesh);
const ReferenceStudy = evolution_mod.reference.ReferenceStudy(SurfaceMesh, ExactInitializer, SurfaceErrorMeasure);

pub const ConfigImpl = struct {
    refinement: u32 = 0,
    steps: u32 = 8,
    dt_scale: f64 = 0.1,
    final_time: f64 = 0.05,
    output_dir: []const u8 = "output/diffusion_surface",
    frames: u32 = 4,

    pub fn timeStep(self: ConfigImpl) f64 {
        std.debug.assert(self.steps > 0);
        return self.final_time / @as(f64, @floatFromInt(self.steps));
    }
};

pub const RunResultImpl = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

pub const ConvergenceResultImpl = struct {
    refinement: u32,
    l2_error: f64,
};

const SurfaceSolveRuntime = struct {
    operator_context: *SurfaceOperatorContext,
    linear_system: linear_system.LinearSystem,
    masses: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const SurfaceMesh,
        dt: f64,
    ) !SurfaceSolveRuntime {
        const operator_context = try SurfaceOperatorContext.init(allocator, mesh);
        errdefer operator_context.deinit();
        const stiffness = (try operator_context.laplacian(0)).stiffness;

        const masses = try assembleLumpedSurfaceMasses(allocator, mesh);
        errdefer allocator.free(masses);

        var assembler = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
        defer assembler.deinit(allocator);

        for (masses, 0..) |mass, row_idx| {
            try assembler.addEntry(allocator, @intCast(row_idx), @intCast(row_idx), mass);

            const row = stiffness.row(@intCast(row_idx));
            for (row.cols, row.vals) |col_idx, value| {
                try assembler.addEntry(allocator, @intCast(row_idx), col_idx, dt * value);
            }
        }

        var system_matrix = try assembler.build(allocator);
        errdefer system_matrix.deinit(allocator);

        var system_runtime = try linear_system.LinearSystem.init(
            allocator,
            system_matrix,
            .{},
        );
        errdefer system_runtime.deinit(allocator);

        return .{
            .operator_context = operator_context,
            .linear_system = system_runtime,
            .masses = masses,
        };
    }

    pub fn deinit(self: *SurfaceSolveRuntime, allocator: std.mem.Allocator) void {
        self.linear_system.deinit(allocator);
        allocator.free(self.masses);
        self.operator_context.deinit();
    }
};

pub const SystemImpl = struct {
    mesh: *const SurfaceMesh,
    diffusion_runtime: SurfaceSolveRuntime,
    temperature: VertexField,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const SurfaceMesh,
        dt: f64,
    ) !SystemImpl {
        var diffusion_runtime = try SurfaceSolveRuntime.init(allocator, mesh, dt);
        errdefer diffusion_runtime.deinit(allocator);

        var temperature = try VertexField.init(allocator, mesh);
        errdefer temperature.deinit(allocator);

        return .{
            .mesh = mesh,
            .diffusion_runtime = diffusion_runtime,
            .temperature = temperature,
        };
    }

    pub fn deinit(self: *SystemImpl, allocator: std.mem.Allocator) void {
        self.temperature.deinit(allocator);
        self.diffusion_runtime.deinit(allocator);
    }

    pub fn values(self: *const SystemImpl) []const f64 {
        return self.temperature.values;
    }
};

pub fn runImpl(
    allocator: std.mem.Allocator,
    config: ConfigImpl,
    writer: anytype,
) !RunResultImpl {
    const result = try simulateCase(allocator, config, writer);
    try writer.print(
        "diffusion_surface: refinement={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.refinement, result.steps, config.timeStep(), result.l2_error },
    );
    return result;
}

pub fn runConvergenceStudyImpl(
    allocator: std.mem.Allocator,
    refinements: []const u32,
) ![]ConvergenceResultImpl {
    const Runner = struct {
        pub fn run(allocator_inner: std.mem.Allocator, refinement: u32) !ConvergenceResultImpl {
            const config = convergenceConfig(refinement);
            const run_result = try simulateCase(allocator_inner, config, null);
            return .{
                .refinement = refinement,
                .l2_error = run_result.l2_error,
            };
        }
    };

    return common.runConvergenceStudy(ConvergenceResultImpl, u32, allocator, refinements, Runner{});
}

fn simulateCase(
    allocator: std.mem.Allocator,
    config: ConfigImpl,
    progress_writer: ?*std.Io.Writer,
) !RunResultImpl {
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.dt_scale > 0.0);
    std.debug.assert(config.final_time > 0.0);

    var mesh = try SurfaceMesh.sphere(allocator, 1.0, config.refinement);
    defer mesh.deinit(allocator);

    const dt = config.timeStep();
    var system = try SystemImpl.init(allocator, &mesh, dt);
    defer system.deinit(allocator);
    initializeState(&mesh, system.temperature.values, 0.0);

    var study = try ReferenceStudy.init(
        allocator,
        &mesh,
        system.values(),
        ExactInitializer{},
        SurfaceErrorMeasure{},
    );
    defer study.deinit(allocator);

    var evolution = try evolution_mod.Evolution(*SystemImpl, SurfaceBackwardEulerMethod).config()
        .dt(dt)
        .init(allocator, &system);
    defer evolution.deinit();

    const loop_result = try common.runReferenceEvolutionLoop(
        SurfaceMesh,
        allocator,
        &mesh,
        &evolution,
        &study,
        .{
            .steps = config.steps,
            .final_time = config.final_time,
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = "diffusion_surface",
            .progress_writer = progress_writer,
        },
    );

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .steps = config.steps,
        .snapshot_count = loop_result.snapshot_count,
        .l2_error = study.l2Error(),
    };
}

fn convergenceConfig(refinement: u32) ConfigImpl {
    const steps = 24 * std.math.pow(u32, 4, refinement);
    return .{
        .refinement = refinement,
        .steps = steps,
        .dt_scale = 0.1,
        .final_time = 0.05,
        .frames = 0,
        .output_dir = "output/diffusion_surface_convergence",
    };
}

fn stepBackwardEuler(
    system: *SystemImpl,
) !void {
    const rhs = system.diffusion_runtime.linear_system.rhsValues();
    const solution = system.diffusion_runtime.linear_system.solutionValues();
    std.debug.assert(rhs.len == system.temperature.values.len);
    std.debug.assert(solution.len == system.temperature.values.len);

    for (rhs, system.temperature.values, system.diffusion_runtime.masses) |*rhs_value, state_value, mass| {
        rhs_value.* = mass * state_value;
    }

    _ = try system.diffusion_runtime.linear_system.solve();

    @memcpy(system.temperature.values, solution);
}

const SurfaceBackwardEulerMethod = struct {
    pub fn initialize(_: std.mem.Allocator, system: *SystemImpl, _: f64) void {
        system.diffusion_runtime.linear_system.seedSolution(system.temperature.values);
    }

    pub fn advance(
        allocator: std.mem.Allocator,
        system: *SystemImpl,
        _: f64,
    ) !void {
        _ = allocator;
        try stepBackwardEuler(system);
    }
};

const ExactInitializer = struct {
    pub fn fill(_: @This(), mesh: *const SurfaceMesh, values: []f64, time: f64) void {
        initializeState(mesh, values, time);
    }
};

const SurfaceErrorMeasure = struct {
    pub fn compute(_: @This(), mesh: *const SurfaceMesh, approx: []const f64, exact: []const f64) f64 {
        return weightedL2Error(mesh, approx, exact);
    }
};
fn assembleLumpedSurfaceMasses(
    allocator: std.mem.Allocator,
    mesh: *const SurfaceMesh,
) ![]f64 {
    const masses = try allocator.alloc(f64, mesh.num_vertices());
    @memset(masses, 0.0);

    const face_vertices = mesh.simplices(2).items(.vertices);
    const face_areas = mesh.simplices(2).items(.volume);
    for (face_vertices, face_areas) |face, area| {
        const lumped = area / 3.0;
        masses[face[0]] += lumped;
        masses[face[1]] += lumped;
        masses[face[2]] += lumped;
    }
    return masses;
}

fn initializeState(
    mesh: *const SurfaceMesh,
    values: []f64,
    time: f64,
) void {
    const coords = mesh.vertices.slice().items(.coords);
    std.debug.assert(values.len == coords.len);
    for (coords, values) |coord, *value| {
        value.* = exactSolution(coord, time);
    }
}

fn exactSolution(coords: [3]f64, time: f64) f64 {
    const eigenvalue = 2.0;
    return std.math.exp(-eigenvalue * time) * coords[2];
}

fn weightedL2Error(mesh: *const SurfaceMesh, approx: []const f64, exact: []const f64) f64 {
    std.debug.assert(approx.len == exact.len);
    const face_vertices = mesh.simplices(2).items(.vertices);
    const face_areas = mesh.simplices(2).items(.volume);
    var error_sq: f64 = 0.0;
    var measure: f64 = 0.0;
    for (face_vertices, face_areas) |face, area| {
        const lumped = area / 3.0;
        const diff0 = approx[face[0]] - exact[face[0]];
        const diff1 = approx[face[1]] - exact[face[1]];
        const diff2 = approx[face[2]] - exact[face[2]];
        error_sq += lumped * (diff0 * diff0 + diff1 * diff1 + diff2 * diff2);
        measure += 3.0 * lumped;
    }
    return @sqrt(error_sq / measure);
}
