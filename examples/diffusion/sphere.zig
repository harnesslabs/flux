const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const cg = flux.math.cg;
const operator_context_mod = flux.operators.context;
const evolution_mod = flux.integrators.evolution;

pub const SurfaceMesh = flux.topology.Mesh(3, 2);
pub const VertexField = flux.forms.Cochain(SurfaceMesh, 0, flux.forms.Primal);
const SurfaceOperatorContext = operator_context_mod.OperatorContext(SurfaceMesh);

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

pub const SurfaceSystem = struct {
    operator_context: *SurfaceOperatorContext,
    system_matrix: sparse.CsrMatrix(f64),
    masses: []f64,
    diagonal: []f64,
    scratch: cg.Scratch,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const SurfaceMesh,
        dt: f64,
    ) !SurfaceSystem {
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

        const diagonal = try allocator.alloc(f64, masses.len);
        errdefer allocator.free(diagonal);
        for (0..masses.len) |row_idx| {
            diagonal[row_idx] = diagonalEntry(system_matrix, @intCast(row_idx));
            std.debug.assert(diagonal[row_idx] > 0.0);
        }

        var scratch = try cg.Scratch.init(allocator, mesh.num_vertices());
        errdefer scratch.deinit(allocator);

        return .{
            .operator_context = operator_context,
            .system_matrix = system_matrix,
            .masses = masses,
            .diagonal = diagonal,
            .scratch = scratch,
        };
    }

    pub fn deinit(self: *SurfaceSystem, allocator: std.mem.Allocator) void {
        self.scratch.deinit(allocator);
        allocator.free(self.diagonal);
        allocator.free(self.masses);
        self.system_matrix.deinit(allocator);
        self.operator_context.deinit();
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
    var system = try SurfaceSystem.init(allocator, &mesh, dt);
    defer system.deinit(allocator);

    var state = try VertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeState(&mesh, state.values, 0.0);

    const Evolution = evolution_mod.ExactEvolution(
        SurfaceMesh,
        f64,
        SurfaceStepperBuilder,
        ExactInitializer,
        SurfaceErrorMeasure,
    );
    var evolution = try Evolution.init(
        allocator,
        &mesh,
        state.values,
        SurfaceStepperBuilder{ .system = &system },
        ExactInitializer{},
        SurfaceErrorMeasure{},
    );
    defer evolution.deinit();

    const loop_result = try common.runExactEvolutionLoop(
        SurfaceMesh,
        allocator,
        &mesh,
        &evolution,
        .{
            .steps = config.steps,
            .dt = dt,
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
        .l2_error = evolution.l2Error(),
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

fn diagonalEntry(matrix: sparse.CsrMatrix(f64), row_idx: u32) f64 {
    const row = matrix.row(row_idx);
    for (row.cols, row.vals) |col_idx, value| {
        if (col_idx == row_idx) return value;
    }
    unreachable;
}

fn stepBackwardEuler(
    system: *SurfaceSystem,
    state_values: []f64,
    rhs: []f64,
    solution: []f64,
) !void {
    std.debug.assert(rhs.len == state_values.len);
    std.debug.assert(solution.len == state_values.len);

    for (rhs, state_values, system.masses) |*rhs_value, state_value, mass| {
        rhs_value.* = mass * state_value;
    }

    var preconditioner = cg.DiagonalPreconditioner{ .diagonal = system.diagonal };
    const solve_result = cg.solve(
        system.system_matrix,
        rhs,
        solution,
        1e-10,
        4000,
        &preconditioner,
        system.scratch,
    );
    if (!solve_result.converged) return error.ConjugateGradientDidNotConverge;

    @memcpy(state_values, solution);
}

const SurfaceStepper = struct {
    system: *SurfaceSystem,
    state_values: []f64,
    rhs: []f64,
    solution: []f64,

    pub fn step(self: *@This(), allocator: std.mem.Allocator) !void {
        _ = allocator;
        try stepBackwardEuler(self.system, self.state_values, self.rhs, self.solution);
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.solution);
        allocator.free(self.rhs);
    }
};

const SurfaceStepperBuilder = struct {
    pub const Stepper = SurfaceStepper;

    system: *SurfaceSystem,

    pub fn initStepper(self: @This(), allocator: std.mem.Allocator, state_values: []f64) !Stepper {
        const len = self.system.masses.len;
        const stepper = Stepper{
            .system = self.system,
            .state_values = state_values,
            .rhs = try allocator.alloc(f64, len),
            .solution = try allocator.alloc(f64, len),
        };
        @memcpy(stepper.solution, state_values);
        return stepper;
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
