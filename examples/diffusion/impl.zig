const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const cg = flux.math.cg;
const exterior_derivative = flux.operators.exterior_derivative;
const hodge_star = flux.operators.hodge_star;
const operator_context_mod = flux.operators.context;

pub const SurfaceKind = enum {
    plane,
    sphere,
};

const SurfaceMesh = flux.topology.Mesh(3, 2);
const SphereVertexField = flux.forms.Cochain(SurfaceMesh, 0, flux.forms.Primal);
const SphereEdgeField = flux.forms.Cochain(SurfaceMesh, 1, flux.forms.Primal);

const SphereConfigImpl = struct {
    refinement: u32 = 0,
    steps: u32 = 8,
    dt_scale: f64 = 0.1,
    final_time: f64 = 0.05,
    output_dir: []const u8 = "output/diffusion_surface",
    frames: u32 = 4,

    pub fn timeStep(self: SphereConfigImpl) f64 {
        std.debug.assert(self.steps > 0);
        return self.final_time / @as(f64, @floatFromInt(self.steps));
    }
};

const SphereRunResultImpl = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

const SphereConvergenceResultImpl = struct {
    refinement: u32,
    l2_error: f64,
};

const SphereSurfaceSystem = struct {
    system_matrix: sparse.CsrMatrix(f64),
    masses: []f64,
    diagonal: []f64,
    scratch: cg.Scratch,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const SurfaceMesh,
        dt: f64,
    ) !SphereSurfaceSystem {
        var stiffness = try assembleSurfaceStiffness(allocator, mesh);
        defer stiffness.deinit(allocator);

        const masses = try assembleLumpedSurfaceMasses(allocator, mesh);
        errdefer allocator.free(masses);

        var assembler = sparse.TripletAssembler(f64).init(
            mesh.num_vertices(),
            mesh.num_vertices(),
        );
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
            diagonal[row_idx] = sphereDiagonalEntry(system_matrix, @intCast(row_idx));
            std.debug.assert(diagonal[row_idx] > 0.0);
        }

        var scratch = try cg.Scratch.init(allocator, mesh.num_vertices());
        errdefer scratch.deinit(allocator);

        return .{
            .system_matrix = system_matrix,
            .masses = masses,
            .diagonal = diagonal,
            .scratch = scratch,
        };
    }

    pub fn deinit(self: *SphereSurfaceSystem, allocator: std.mem.Allocator) void {
        self.scratch.deinit(allocator);
        allocator.free(self.diagonal);
        allocator.free(self.masses);
        self.system_matrix.deinit(allocator);
    }
};

fn runSphereImpl(
    allocator: std.mem.Allocator,
    config: SphereConfigImpl,
    writer: anytype,
) !SphereRunResultImpl {
    const result = try sphereSimulateCase(allocator, config);
    try writer.print(
        "diffusion_surface: refinement={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.refinement, result.steps, config.timeStep(), result.l2_error },
    );
    return result;
}

fn runSphereConvergenceStudyImpl(
    allocator: std.mem.Allocator,
    refinements: []const u32,
) ![]SphereConvergenceResultImpl {
    const Runner = struct {
        pub fn run(allocator_inner: std.mem.Allocator, refinement: u32) !SphereConvergenceResultImpl {
            const config = sphereConvergenceConfig(refinement);
            const run_result = try sphereSimulateCase(allocator_inner, config);
            return .{
                .refinement = refinement,
                .l2_error = run_result.l2_error,
            };
        }
    };

    return common.runConvergenceStudy(SphereConvergenceResultImpl, u32, allocator, refinements, Runner{});
}

fn sphereSimulateCase(
    allocator: std.mem.Allocator,
    config: SphereConfigImpl,
) !SphereRunResultImpl {
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.dt_scale > 0.0);
    std.debug.assert(config.final_time > 0.0);

    var mesh = try buildSphereMesh(allocator, config.refinement);
    defer mesh.deinit(allocator);

    const dt = config.timeStep();
    var system = try SphereSurfaceSystem.init(allocator, &mesh, dt);
    defer system.deinit(allocator);

    var state = try SphereVertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeSphereState(&mesh, state.values, 0.0);

    const exact = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(exact);

    const rhs = try allocator.alloc(f64, state.values.len);
    defer allocator.free(rhs);
    const solution = try allocator.alloc(f64, state.values.len);
    defer allocator.free(solution);
    @memcpy(solution, state.values);

    const loop_result = try common.runExactFieldLoop(
        SurfaceMesh,
        allocator,
        &mesh,
        state.values,
        exact,
        .{
            .steps = config.steps,
            .dt = dt,
            .final_time = config.final_time,
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = "diffusion_surface",
        },
        SphereSurfaceExactInitializer{},
        SphereSurfaceStepper{
            .system = &system,
            .state_values = state.values,
            .rhs = rhs,
            .solution = solution,
        },
    );

    const l2_error = sphereWeightedL2Error(&mesh, state.values, exact);

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .steps = config.steps,
        .snapshot_count = loop_result.snapshot_count,
        .l2_error = l2_error,
    };
}

fn sphereConvergenceConfig(refinement: u32) SphereConfigImpl {
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

fn buildSphereMesh(
    allocator: std.mem.Allocator,
    refinement: u32,
) !SurfaceMesh {
    var polyhedron = try buildRefinedOctahedron(allocator, refinement);
    defer polyhedron.deinit(allocator);

    const oriented_faces = try allocator.dupe([3]u32, polyhedron.faces);
    defer allocator.free(oriented_faces);
    orientFacesOutward(polyhedron.vertices, oriented_faces);

    return SurfaceMesh.from_triangles(allocator, polyhedron.vertices, oriented_faces);
}

const Polyhedron = struct {
    vertices: [][3]f64,
    faces: [][3]u32,

    pub fn deinit(self: *Polyhedron, allocator: std.mem.Allocator) void {
        allocator.free(self.faces);
        allocator.free(self.vertices);
    }
};

fn buildRefinedOctahedron(
    allocator: std.mem.Allocator,
    refinement: u32,
) !Polyhedron {
    var polyhedron = Polyhedron{
        .vertices = try allocator.dupe([3]f64, &initial_octahedron_vertices),
        .faces = try allocator.dupe([3]u32, &initial_octahedron_faces),
    };
    errdefer polyhedron.deinit(allocator);

    for (0..refinement) |_| {
        const next = try refinePolyhedron(allocator, polyhedron);
        polyhedron.deinit(allocator);
        polyhedron = next;
    }

    return polyhedron;
}

fn refinePolyhedron(
    allocator: std.mem.Allocator,
    polyhedron: Polyhedron,
) !Polyhedron {
    var vertices = try std.ArrayList([3]f64).initCapacity(allocator, polyhedron.vertices.len + polyhedron.faces.len);
    defer vertices.deinit(allocator);
    try vertices.appendSlice(allocator, polyhedron.vertices);

    var faces = try std.ArrayList([3]u32).initCapacity(allocator, polyhedron.faces.len * 4);
    defer faces.deinit(allocator);

    var midpoint_map = std.AutoHashMap([2]u32, u32).init(allocator);
    defer midpoint_map.deinit();

    for (polyhedron.faces) |face| {
        const ab = try midpointIndex(allocator, &vertices, &midpoint_map, face[0], face[1]);
        const bc = try midpointIndex(allocator, &vertices, &midpoint_map, face[1], face[2]);
        const ca = try midpointIndex(allocator, &vertices, &midpoint_map, face[2], face[0]);

        try faces.append(allocator, .{ face[0], ab, ca });
        try faces.append(allocator, .{ face[1], bc, ab });
        try faces.append(allocator, .{ face[2], ca, bc });
        try faces.append(allocator, .{ ab, bc, ca });
    }

    return .{
        .vertices = try vertices.toOwnedSlice(allocator),
        .faces = try faces.toOwnedSlice(allocator),
    };
}

fn midpointIndex(
    allocator: std.mem.Allocator,
    vertices: *std.ArrayList([3]f64),
    midpoint_map: *std.AutoHashMap([2]u32, u32),
    a: u32,
    b: u32,
) !u32 {
    const key = canonicalEdge(a, b);
    const gop = try midpoint_map.getOrPut(key);
    if (gop.found_existing) return gop.value_ptr.*;

    const midpoint = normalize3(add3(vertices.items[a], vertices.items[b]));
    const new_index: u32 = @intCast(vertices.items.len);
    try vertices.append(allocator, midpoint);
    gop.value_ptr.* = new_index;
    return new_index;
}

fn orientFacesOutward(
    embedded_vertices: []const [3]f64,
    faces: [][3]u32,
) void {
    for (faces) |*face| {
        const a = embedded_vertices[face.*[0]];
        const b = embedded_vertices[face.*[1]];
        const c = embedded_vertices[face.*[2]];
        const ab = sub3(b, a);
        const ac = sub3(c, a);
        const normal = cross3(ab, ac);
        const centroid = normalize3(.{
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        });
        if (dot3(normal, centroid) < 0.0) {
            const tmp = face.*[1];
            face.*[1] = face.*[2];
            face.*[2] = tmp;
        }
    }
}

fn assembleSurfaceStiffness(
    allocator: std.mem.Allocator,
    mesh: *const SurfaceMesh,
) !sparse.CsrMatrix(f64) {
    var assembler = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
    defer assembler.deinit(allocator);

    var basis = try SphereVertexField.init(allocator, mesh);
    defer basis.deinit(allocator);

    const boundary_1 = mesh.boundary(1);
    const column = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(column);

    for (0..mesh.num_vertices()) |vertex_idx| {
        @memset(basis.values, 0.0);
        basis.values[vertex_idx] = 1.0;

        var gradient = try exterior_derivative.exterior_derivative(allocator, basis);
        defer gradient.deinit(allocator);

        var starred = try hodge_star.hodge_star(allocator, gradient);
        defer starred.deinit(allocator);

        @memset(column, 0.0);
        boundary_1.transpose_multiply(starred.values, column);
        for (column, 0..) |value, row_idx| {
            if (@abs(value) < 1e-13) continue;
            try assembler.addEntry(allocator, @intCast(row_idx), @intCast(vertex_idx), value);
        }
    }

    return assembler.build(allocator);
}

fn sphereDiagonalEntry(matrix: sparse.CsrMatrix(f64), row_idx: u32) f64 {
    const row = matrix.row(row_idx);
    for (row.cols, row.vals) |col_idx, value| {
        if (col_idx == row_idx) return value;
    }
    unreachable;
}

fn sphereStepBackwardEuler(
    system: *SphereSurfaceSystem,
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

const SphereSurfaceStepper = struct {
    system: *SphereSurfaceSystem,
    state_values: []f64,
    rhs: []f64,
    solution: []f64,

    pub fn step(self: @This(), allocator: std.mem.Allocator) !void {
        _ = allocator;
        try sphereStepBackwardEuler(self.system, self.state_values, self.rhs, self.solution);
    }
};

const SphereSurfaceExactInitializer = struct {
    pub fn fill(_: @This(), mesh: *const SurfaceMesh, values: []f64, time: f64) void {
        initializeSphereState(mesh, values, time);
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

fn initializeSphereState(
    mesh: *const SurfaceMesh,
    values: []f64,
    time: f64,
) void {
    const coords = mesh.vertices.slice().items(.coords);
    std.debug.assert(values.len == coords.len);
    for (coords, values) |coord, *value| {
        value.* = sphereExactSolution(coord, time);
    }
}

fn sphereExactSolution(coords: [3]f64, time: f64) f64 {
    const eigenvalue = 2.0;
    return std.math.exp(-eigenvalue * time) * coords[2];
}

fn sphereWeightedL2Error(mesh: *const SurfaceMesh, approx: []const f64, exact: []const f64) f64 {
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

fn canonicalEdge(a: u32, b: u32) [2]u32 {
    return if (a < b) .{ a, b } else .{ b, a };
}

fn add3(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

fn sub3(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

fn dot3(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

fn norm3(a: [3]f64) f64 {
    return @sqrt(dot3(a, a));
}

fn normalize3(a: [3]f64) [3]f64 {
    const inv_norm = 1.0 / norm3(a);
    return .{ a[0] * inv_norm, a[1] * inv_norm, a[2] * inv_norm };
}

fn cross3(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

const initial_octahedron_vertices = [_][3]f64{
    .{ 1.0, 0.0, 0.0 },
    .{ -1.0, 0.0, 0.0 },
    .{ 0.0, 1.0, 0.0 },
    .{ 0.0, -1.0, 0.0 },
    .{ 0.0, 0.0, 1.0 },
    .{ 0.0, 0.0, -1.0 },
};

const initial_octahedron_faces = [_][3]u32{
    .{ 0, 4, 2 },
    .{ 2, 4, 1 },
    .{ 1, 4, 3 },
    .{ 3, 4, 0 },
    .{ 2, 5, 0 },
    .{ 1, 5, 2 },
    .{ 3, 5, 1 },
    .{ 0, 5, 3 },
};

test "surface diffusion error decreases under sphere refinement" {
    const allocator = testing.allocator;
    const refinements = [_]u32{ 1, 2, 3 };

    const results = try runSphereConvergenceStudyImpl(allocator, &refinements);
    defer allocator.free(results);

    try testing.expectEqual(refinements.len, results.len);
    try testing.expect(results[1].l2_error < results[0].l2_error);
    try testing.expect(results[2].l2_error < results[1].l2_error);
    try testing.expect(results[2].l2_error < 1e-3);
}

test "surface diffusion example runs as standalone binary configuration" {
    const allocator = testing.allocator;
    var log_buffer = std.ArrayListUnmanaged(u8){};
    defer log_buffer.deinit(allocator);

    const result = try runSphereImpl(allocator, .{
        .refinement = 0,
        .steps = 2,
        .frames = 0,
        .output_dir = "output/diffusion-surface-test",
    }, log_buffer.writer(allocator));

    try testing.expectEqual(@as(u32, 2), result.steps);
    try testing.expect(result.elapsed_s >= 0.0);
}

const Mesh2D = flux.topology.Mesh(2, 2);
const VertexField = flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal);

const convergence_time = 0.02;

const InitialCondition = enum {
    zero,
    sine_mode,
};

const PlaneConfigImpl = struct {
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

    pub fn dt(self: PlaneConfigImpl) f64 {
        if (self.dt_override) |value| return value;
        const h = self.domain / @as(f64, @floatFromInt(self.grid));
        return self.dt_scale * h * h;
    }
};

const PlaneRunResultImpl = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

const PlaneConvergenceResultImpl = struct {
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

fn runPlaneImpl(
    allocator: std.mem.Allocator,
    config: PlaneConfigImpl,
    writer: anytype,
) !PlaneRunResultImpl {
    const result = try simulateCase(allocator, config, .sine_mode, null);
    try writer.print(
        "heat: grid={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.grid, config.steps, config.dt(), result.l2_error },
    );
    return result;
}

fn runPlaneConvergenceStudyImpl(
    allocator: std.mem.Allocator,
    grids: []const u32,
) ![]PlaneConvergenceResultImpl {
    const Runner = struct {
        pub fn run(allocator_inner: std.mem.Allocator, grid: u32) !PlaneConvergenceResultImpl {
            const config = convergenceConfig(grid);
            const run_result = try simulateCase(allocator_inner, config, .sine_mode, convergence_time);
            return .{
                .grid = grid,
                .l2_error = run_result.l2_error,
            };
        }
    };

    return common.runConvergenceStudy(PlaneConvergenceResultImpl, u32, allocator, grids, Runner{});
}

fn simulateCase(
    allocator: std.mem.Allocator,
    config: PlaneConfigImpl,
    initial_condition: InitialCondition,
    final_time_override: ?f64,
) !PlaneRunResultImpl {
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

fn convergenceConfig(grid: u32) PlaneConfigImpl {
    const probe = PlaneConfigImpl{ .grid = grid };
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
    const config = PlaneConfigImpl{
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
    const config = PlaneConfigImpl{
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

    const result = try run(.plane, allocator, .{
        .grid = 4,
        .steps = 2,
        .frames = 0,
        .output_dir = "output/heat-test",
    }, log_buffer.writer(allocator));

    try testing.expectEqual(@as(u32, 2), result.steps);
    try testing.expect(result.elapsed_s >= 0.0);
}

pub fn Mesh(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => Mesh2D,
        .sphere => SurfaceMesh,
    };
}

pub fn State(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => VertexField,
        .sphere => SphereVertexField,
    };
}

pub fn Config(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => PlaneConfigImpl,
        .sphere => SphereConfigImpl,
    };
}

pub fn RunResult(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => PlaneRunResultImpl,
        .sphere => SphereRunResultImpl,
    };
}

pub fn ConvergenceResult(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => PlaneConvergenceResultImpl,
        .sphere => SphereConvergenceResultImpl,
    };
}

pub fn run(comptime surface_kind: SurfaceKind, allocator: std.mem.Allocator, config: Config(surface_kind), writer: anytype) !RunResult(surface_kind) {
    return switch (surface_kind) {
        .plane => try runPlaneImpl(allocator, config, writer),
        .sphere => try runSphereImpl(allocator, config, writer),
    };
}

pub fn runConvergenceStudy(comptime surface_kind: SurfaceKind, allocator: std.mem.Allocator, params: []const u32) ![]ConvergenceResult(surface_kind) {
    return switch (surface_kind) {
        .plane => try runPlaneConvergenceStudyImpl(allocator, params),
        .sphere => try runSphereConvergenceStudyImpl(allocator, params),
    };
}
