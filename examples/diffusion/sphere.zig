const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const sparse = flux.math.sparse;
const cg = flux.math.cg;
const exterior_derivative = flux.operators.exterior_derivative;
const hodge_star = flux.operators.hodge_star;

pub const SurfaceMesh = flux.topology.Mesh(3, 2);
pub const VertexField = flux.forms.Cochain(SurfaceMesh, 0, flux.forms.Primal);

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

const SurfaceSystem = struct {
    system_matrix: sparse.CsrMatrix(f64),
    masses: []f64,
    diagonal: []f64,
    scratch: cg.Scratch,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const SurfaceMesh,
        dt: f64,
    ) !SurfaceSystem {
        var stiffness = try assembleSurfaceStiffness(allocator, mesh);
        defer stiffness.deinit(allocator);

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

    var mesh = try buildSphereMesh(allocator, config.refinement);
    defer mesh.deinit(allocator);

    const dt = config.timeStep();
    var system = try SurfaceSystem.init(allocator, &mesh, dt);
    defer system.deinit(allocator);

    var state = try VertexField.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeState(&mesh, state.values, 0.0);

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
            .progress_writer = progress_writer,
        },
        ExactInitializer{},
        Stepper{
            .system = &system,
            .state_values = state.values,
            .rhs = rhs,
            .solution = solution,
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

    var basis = try VertexField.init(allocator, mesh);
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

const Stepper = struct {
    system: *SurfaceSystem,
    state_values: []f64,
    rhs: []f64,
    solution: []f64,

    pub fn step(self: @This(), allocator: std.mem.Allocator) !void {
        _ = allocator;
        try stepBackwardEuler(self.system, self.state_values, self.rhs, self.solution);
    }
};

const ExactInitializer = struct {
    pub fn fill(_: @This(), mesh: *const SurfaceMesh, values: []f64, time: f64) void {
        initializeState(mesh, values, time);
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
