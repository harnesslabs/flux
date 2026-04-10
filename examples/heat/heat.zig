const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const common = @import("examples_common");

pub const Mesh2D = flux.topology.Mesh(2, 2);
pub const VertexField = flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal);
const sparse = flux.math.sparse;
const cg = flux.math.cg;
const operator_context_mod = flux.operators.context;

const convergence_time = 0.02;

const InitialCondition = enum {
    zero,
    sine_mode,
};

pub const Config = struct {
    grid: u32 = 8,
    steps: u32 = 8,
    domain: f64 = 1.0,
    dt_scale: f64 = 0.1,
    /// Optional explicit timestep override. When set, takes precedence over
    /// the `dt_scale * h^2` parabolic default. Plumbed through from the
    /// shared `--dt` CLI flag.
    dt_override: ?f64 = null,
    output_dir: []const u8 = "output/heat",
    frames: u32 = 4,

    pub fn dt(self: Config) f64 {
        if (self.dt_override) |value| return value;
        const h = self.domain / @as(f64, @floatFromInt(self.grid));
        return self.dt_scale * h * h;
    }
};

pub const RunResult = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

pub const ConvergenceResult = struct {
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

pub fn run(
    allocator: std.mem.Allocator,
    config: Config,
    writer: anytype,
) !RunResult {
    const result = try simulateCase(allocator, config, .sine_mode, null);
    try writer.print(
        "heat: grid={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.grid, config.steps, config.dt(), result.l2_error },
    );
    return result;
}

pub fn runConvergenceStudy(
    allocator: std.mem.Allocator,
    grids: []const u32,
) ![]ConvergenceResult {
    const results = try allocator.alloc(ConvergenceResult, grids.len);
    errdefer allocator.free(results);

    for (grids, 0..) |grid, idx| {
        const config = convergenceConfig(grid);
        const run_result = try simulateCase(allocator, config, .sine_mode, convergence_time);
        results[idx] = .{
            .grid = grid,
            .l2_error = run_result.l2_error,
        };
    }

    return results;
}

fn simulateCase(
    allocator: std.mem.Allocator,
    config: Config,
    initial_condition: InitialCondition,
    final_time_override: ?f64,
) !RunResult {
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

    const exact = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(exact);

    const reduced_rhs = try allocator.alloc(f64, heat_system.interior_vertices.len);
    defer allocator.free(reduced_rhs);
    const reduced_solution = try allocator.alloc(f64, heat_system.interior_vertices.len);
    defer allocator.free(reduced_solution);
    seedReducedSolution(&heat_system, state.values, reduced_solution);

    const loop_result = try common.runExactFieldLoop(
        Mesh2D,
        allocator,
        &mesh,
        state.values,
        exact,
        .{
            .steps = config.steps,
            .dt = dt,
            .final_time = total_time,
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = "heat",
        },
        HeatExactInitializer{ .initial_condition = initial_condition },
        HeatStepper{
            .mesh = &mesh,
            .heat_system = &heat_system,
            .state_values = state.values,
            .reduced_rhs = reduced_rhs,
            .reduced_solution = reduced_solution,
        },
    );

    const l2_error = weightedL2Error(&mesh, state.values, exact);

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .steps = config.steps,
        .snapshot_count = loop_result.snapshot_count,
        .l2_error = l2_error,
    };
}

fn convergenceConfig(grid: u32) Config {
    const probe = Config{ .grid = grid };
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

const HeatStepper = struct {
    mesh: *const Mesh2D,
    heat_system: *const HeatSystem,
    state_values: []f64,
    reduced_rhs: []f64,
    reduced_solution: []f64,

    pub fn step(self: @This(), allocator: std.mem.Allocator) !void {
        try stepBackwardEuler(
            allocator,
            self.mesh,
            self.heat_system,
            self.state_values,
            self.reduced_rhs,
            self.reduced_solution,
        );
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
        if (is_boundary) {
            state_values[vertex_idx] = 0.0;
        }
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

test "heat zero initial data stays zero under backward Euler" {
    const allocator = testing.allocator;
    const config = Config{
        .grid = 4,
        .steps = 3,
        .frames = 0,
    };
    const result = try simulateCase(allocator, config, .zero, 0.01);

    try testing.expectApproxEqAbs(@as(f64, 0.0), result.l2_error, 1e-12);
}

test "heat convergence study reaches second-order spatial rate" {
    const allocator = testing.allocator;
    const grids = [_]u32{ 8, 16, 32 };
    const results = try runConvergenceStudy(allocator, &grids);
    defer allocator.free(results);

    try testing.expectEqual(grids.len, results.len);
    for (0..results.len - 1) |idx| {
        const rate = std.math.log(f64, 2.0, results[idx].l2_error / results[idx + 1].l2_error);
        try testing.expect(rate > 1.75);
    }
}

test "heat example preserves homogeneous Dirichlet boundary values" {
    const allocator = testing.allocator;
    const config = Config{
        .grid = 6,
        .steps = 4,
        .frames = 0,
    };

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);
    const operator_context = try operator_context_mod.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();

    const dt = config.dt();
    var heat_system = try HeatSystem.init(allocator, &mesh, operator_context, dt);
    defer heat_system.deinit(allocator);

    var state = try VertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeState(&mesh, state.values, .sine_mode, 0.0);

    const reduced_rhs = try allocator.alloc(f64, heat_system.interior_vertices.len);
    defer allocator.free(reduced_rhs);
    const reduced_solution = try allocator.alloc(f64, heat_system.interior_vertices.len);
    defer allocator.free(reduced_solution);
    seedReducedSolution(&heat_system, state.values, reduced_solution);

    try stepBackwardEuler(allocator, &mesh, &heat_system, state.values, reduced_rhs, reduced_solution);
    for (heat_system.boundary_mask, state.values) |is_boundary, value| {
        if (is_boundary) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-12);
        }
    }
}

test "heat example runs as standalone binary configuration" {
    const allocator = testing.allocator;
    var log_buffer = std.ArrayListUnmanaged(u8){};
    defer log_buffer.deinit(allocator);

    const result = try run(allocator, .{
        .grid = 4,
        .steps = 2,
        .frames = 0,
        .output_dir = "output/heat-test",
    }, log_buffer.writer(allocator));

    try testing.expectEqual(@as(u32, 2), result.steps);
    try testing.expect(result.elapsed_s >= 0.0);
}
