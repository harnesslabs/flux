const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const std_fs = std.fs;

const flux_io = flux.io;
const sparse = flux.math.sparse;
const cg = flux.math.cg;
const exterior_derivative = flux.operators.exterior_derivative;
const hodge_star = flux.operators.hodge_star;

pub const ReferenceMesh = flux.Mesh(2, 2);
pub const EmbeddedMesh = flux.Mesh(3, 2);
pub const VertexField = flux.Cochain(ReferenceMesh, 0, flux.Primal);
pub const EdgeField = flux.Cochain(ReferenceMesh, 1, flux.Primal);
pub const Metric2D = hodge_star.Metric(ReferenceMesh, .riemannian);

pub const Config = struct {
    refinement: u32 = 0,
    steps: u32 = 8,
    dt_scale: f64 = 0.1,
    final_time: f64 = 0.05,
    output_dir: []const u8 = "output/diffusion_surface",
    frames: u32 = 4,

    pub fn timeStep(self: Config) f64 {
        std.debug.assert(self.steps > 0);
        return self.final_time / @as(f64, @floatFromInt(self.steps));
    }
};

pub const RunResult = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

pub const ConvergenceResult = struct {
    refinement: u32,
    l2_error: f64,
};

// TODO(#154): This struct is a workaround. It exists because Mesh(D,K) with
// D > K does not yet honor the metric induced by its own embedding — so we
// must carry a parallel Mesh(2,2) reference mesh, an explicit per-face
// metric tensor array, and the embedded Mesh(3,2) side by side. Once #154
// lands, this collapses to a single `Mesh(3, 2)` and the stereographic
// projection / metric tensor machinery in this file can be deleted. The
// long-term endgame is the truly intrinsic IntrinsicMesh(K) tracked in
// project/horizons.md.
const SphereGeometry = struct {
    reference_mesh: ReferenceMesh,
    embedded_mesh: EmbeddedMesh,
    metric_tensors: [][2][2]f64,

    pub fn deinit(self: *SphereGeometry, allocator: std.mem.Allocator) void {
        allocator.free(self.metric_tensors);
        self.embedded_mesh.deinit(allocator);
        self.reference_mesh.deinit(allocator);
    }
};

const SurfaceSystem = struct {
    system_matrix: sparse.CsrMatrix(f64),
    masses: []f64,
    diagonal: []f64,
    scratch: cg.Scratch,

    pub fn init(
        allocator: std.mem.Allocator,
        geometry: *const SphereGeometry,
        dt: f64,
    ) !SurfaceSystem {
        const metric = Metric2D{
            .top_simplex_tensors = geometry.metric_tensors,
        };

        var stiffness = try assembleMetricStiffness(allocator, &geometry.reference_mesh, metric);
        defer stiffness.deinit(allocator);

        const masses = try assembleLumpedSurfaceMasses(allocator, &geometry.embedded_mesh);
        errdefer allocator.free(masses);

        var assembler = sparse.TripletAssembler(f64).init(
            geometry.reference_mesh.num_vertices(),
            geometry.reference_mesh.num_vertices(),
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
            diagonal[row_idx] = diagonalEntry(system_matrix, @intCast(row_idx));
            std.debug.assert(diagonal[row_idx] > 0.0);
        }

        var scratch = try cg.Scratch.init(allocator, geometry.reference_mesh.num_vertices());
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

pub fn run(
    allocator: std.mem.Allocator,
    config: Config,
    writer: anytype,
) !RunResult {
    const result = try simulateCase(allocator, config);
    try writer.print(
        "diffusion_surface: refinement={d} steps={d} dt={d:.6} l2_error={e}\n",
        .{ config.refinement, result.steps, config.timeStep(), result.l2_error },
    );
    return result;
}

pub fn runConvergenceStudy(
    allocator: std.mem.Allocator,
    refinements: []const u32,
) ![]ConvergenceResult {
    const results = try allocator.alloc(ConvergenceResult, refinements.len);
    errdefer allocator.free(results);

    for (refinements, 0..) |refinement, idx| {
        const config = convergenceConfig(refinement);
        const run_result = try simulateCase(allocator, config);
        results[idx] = .{
            .refinement = refinement,
            .l2_error = run_result.l2_error,
        };
    }

    return results;
}

fn simulateCase(
    allocator: std.mem.Allocator,
    config: Config,
) !RunResult {
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.dt_scale > 0.0);
    std.debug.assert(config.final_time > 0.0);

    var geometry = try buildSphereGeometry(allocator, config.refinement);
    defer geometry.deinit(allocator);

    const dt = config.timeStep();
    var system = try SurfaceSystem.init(allocator, &geometry, dt);
    defer system.deinit(allocator);

    var state = try VertexField.init(allocator, &geometry.reference_mesh);
    defer state.deinit(allocator);
    initializeState(&geometry.embedded_mesh, state.values, 0.0);

    const exact = try allocator.alloc(f64, geometry.reference_mesh.num_vertices());
    defer allocator.free(exact);

    const has_output = config.frames > 0;
    const interval = if (has_output) outputInterval(config) else 1;
    const snapshot_capacity: u32 = if (has_output) (config.steps / interval) + 2 else 0;
    var pvd_entries: []flux_io.PvdEntry = &.{};
    var filename_bufs: [][flux_io.max_snapshot_filename_length]u8 = &.{};
    if (has_output) {
        try ensureDir(config.output_dir);
        pvd_entries = try allocator.alloc(flux_io.PvdEntry, snapshot_capacity);
        filename_bufs = try allocator.alloc([flux_io.max_snapshot_filename_length]u8, snapshot_capacity);
    }
    defer if (has_output) {
        allocator.free(filename_bufs);
        allocator.free(pvd_entries);
    };

    var snapshot_count: u32 = 0;
    if (has_output) {
        initializeState(&geometry.embedded_mesh, exact, 0.0);
        try writeSnapshot(
            allocator,
            config.output_dir,
            snapshotName(&filename_bufs, snapshot_count),
            &geometry.embedded_mesh,
            state.values,
            exact,
        );
        pvd_entries[snapshot_count] = .{
            .timestep = 0.0,
            .filename = snapshotName(&filename_bufs, snapshot_count),
        };
        snapshot_count += 1;
    }

    const rhs = try allocator.alloc(f64, state.values.len);
    defer allocator.free(rhs);
    const solution = try allocator.alloc(f64, state.values.len);
    defer allocator.free(solution);
    @memcpy(solution, state.values);

    const start_ns = std.time.nanoTimestamp();
    for (0..config.steps) |step_idx| {
        const next_time = dt * @as(f64, @floatFromInt(step_idx + 1));
        try stepBackwardEuler(&system, state.values, rhs, solution);
        if (has_output and ((step_idx + 1) % interval == 0 or step_idx + 1 == config.steps)) {
            initializeState(&geometry.embedded_mesh, exact, next_time);
            try writeSnapshot(
                allocator,
                config.output_dir,
                snapshotName(&filename_bufs, snapshot_count),
                &geometry.embedded_mesh,
                state.values,
                exact,
            );
            pvd_entries[snapshot_count] = .{
                .timestep = next_time,
                .filename = snapshotName(&filename_bufs, snapshot_count),
            };
            snapshot_count += 1;
        }
    }
    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

    initializeState(&geometry.embedded_mesh, exact, config.final_time);
    const l2_error = weightedL2Error(&geometry.embedded_mesh, state.values, exact);

    if (has_output and snapshot_count > 0) {
        try writePvd(allocator, config.output_dir, "diffusion_surface", pvd_entries[0..snapshot_count]);
    }

    return .{
        .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
        .steps = config.steps,
        .snapshot_count = snapshot_count,
        .l2_error = l2_error,
    };
}

fn convergenceConfig(refinement: u32) Config {
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

fn outputInterval(config: Config) u32 {
    std.debug.assert(config.frames > 0);
    return @max(1, config.steps / config.frames);
}

fn buildSphereGeometry(
    allocator: std.mem.Allocator,
    refinement: u32,
) !SphereGeometry {
    var polyhedron = try buildRefinedOctahedron(allocator, refinement);
    defer polyhedron.deinit(allocator);

    const projection = projectionFrame(polyhedron.vertices);
    const projected_vertices = try allocator.alloc([2]f64, polyhedron.vertices.len);
    defer allocator.free(projected_vertices);
    for (polyhedron.vertices, projected_vertices) |vertex, *projected| {
        projected.* = stereographicProject(vertex, projection);
    }

    const oriented_faces = try allocator.dupe([3]u32, polyhedron.faces);
    defer allocator.free(oriented_faces);
    orientFaces(projected_vertices, oriented_faces);

    var reference_mesh = try ReferenceMesh.from_triangles(allocator, projected_vertices, oriented_faces);
    errdefer reference_mesh.deinit(allocator);

    var embedded_mesh = try EmbeddedMesh.from_triangles(allocator, polyhedron.vertices, oriented_faces);
    errdefer embedded_mesh.deinit(allocator);

    const metric_tensors = try allocator.alloc([2][2]f64, oriented_faces.len);
    errdefer allocator.free(metric_tensors);
    for (oriented_faces, 0..) |face, face_idx| {
        metric_tensors[face_idx] = metricTensorForFace(projected_vertices, polyhedron.vertices, face);
    }

    return .{
        .reference_mesh = reference_mesh,
        .embedded_mesh = embedded_mesh,
        .metric_tensors = metric_tensors,
    };
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

const ProjectionFrame = struct {
    pole: [3]f64,
    tangent_x: [3]f64,
    tangent_y: [3]f64,
};

fn projectionFrame(vertices: []const [3]f64) ProjectionFrame {
    var best_pole = normalize3(.{ 1.0, 2.0, 3.0 });
    var best_margin: f64 = -std.math.inf(f64);

    for (projection_candidates) |candidate| {
        const pole = normalize3(candidate);
        var margin = std.math.inf(f64);
        for (vertices) |vertex| {
            margin = @min(margin, 1.0 - dot3(vertex, pole));
        }
        if (margin > best_margin) {
            best_margin = margin;
            best_pole = pole;
        }
    }

    const reference = if (@abs(best_pole[2]) < 0.9) [3]f64{ 0.0, 0.0, 1.0 } else [3]f64{ 1.0, 0.0, 0.0 };
    const tangent_x = normalize3(cross3(reference, best_pole));
    const tangent_y = cross3(best_pole, tangent_x);
    return .{
        .pole = best_pole,
        .tangent_x = tangent_x,
        .tangent_y = tangent_y,
    };
}

fn stereographicProject(vertex: [3]f64, frame: ProjectionFrame) [2]f64 {
    const denom = 1.0 - dot3(vertex, frame.pole);
    std.debug.assert(denom > 1e-6);
    return .{
        dot3(vertex, frame.tangent_x) / denom,
        dot3(vertex, frame.tangent_y) / denom,
    };
}

fn orientFaces(
    reference_vertices: []const [2]f64,
    faces: [][3]u32,
) void {
    for (faces) |*face| {
        const area = signedTriangleArea(
            reference_vertices[face.*[0]],
            reference_vertices[face.*[1]],
            reference_vertices[face.*[2]],
        );
        if (area < 0.0) {
            const tmp = face.*[1];
            face.*[1] = face.*[2];
            face.*[2] = tmp;
        }
    }
}

fn metricTensorForFace(
    reference_vertices: []const [2]f64,
    embedded_vertices: []const [3]f64,
    face: [3]u32,
) [2][2]f64 {
    const x0 = reference_vertices[face[0]];
    const x1 = reference_vertices[face[1]];
    const x2 = reference_vertices[face[2]];
    const y0 = embedded_vertices[face[0]];
    const y1 = embedded_vertices[face[1]];
    const y2 = embedded_vertices[face[2]];

    const b11 = x1[0] - x0[0];
    const b12 = x2[0] - x0[0];
    const b21 = x1[1] - x0[1];
    const b22 = x2[1] - x0[1];
    const det_b = b11 * b22 - b12 * b21;
    std.debug.assert(@abs(det_b) > 1e-12);

    const inv_b = [2][2]f64{
        .{ b22 / det_b, -b12 / det_b },
        .{ -b21 / det_b, b11 / det_b },
    };

    const v1 = sub3(y1, y0);
    const v2 = sub3(y2, y0);
    const gram = [2][2]f64{
        .{ dot3(v1, v1), dot3(v1, v2) },
        .{ dot3(v2, v1), dot3(v2, v2) },
    };

    return mulMat2(transpose2(inv_b), mulMat2(gram, inv_b));
}

fn assembleMetricStiffness(
    allocator: std.mem.Allocator,
    mesh: *const ReferenceMesh,
    metric: Metric2D,
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

        var starred = try hodge_star.hodge_star_with_metric(allocator, metric, gradient);
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

fn assembleLumpedSurfaceMasses(
    allocator: std.mem.Allocator,
    mesh: *const EmbeddedMesh,
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
    mesh: *const EmbeddedMesh,
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

fn weightedL2Error(mesh: *const EmbeddedMesh, approx: []const f64, exact: []const f64) f64 {
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

fn snapshotName(
    filename_bufs: *const [][flux_io.max_snapshot_filename_length]u8,
    snapshot_idx: u32,
) []const u8 {
    return flux_io.snapshot_filename(
        @constCast(&filename_bufs.*[snapshot_idx]),
        "diffusion_surface",
        snapshot_idx,
    );
}

fn writeSnapshot(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    filename: []const u8,
    mesh: *const EmbeddedMesh,
    state: []const f64,
    exact: []const f64,
) !void {
    const error_values = try allocator.alloc(f64, state.len);
    defer allocator.free(error_values);
    for (state, exact, error_values) |approx_value, exact_value, *error_value| {
        error_value.* = approx_value - exact_value;
    }

    const point_data = [_]flux_io.DataArraySlice{
        .{ .name = "temperature", .values = state },
        .{ .name = "temperature_exact", .values = exact },
        .{ .name = "temperature_error", .values = error_values },
    };

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);
    try flux_io.write(output.writer(allocator), 3, 2, mesh.*, &point_data, &.{});

    var dir = try std_fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(filename, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn writePvd(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    base_name: []const u8,
    entries: []const flux_io.PvdEntry,
) !void {
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);
    try flux_io.write_pvd(output.writer(allocator), entries);

    var pvd_buf: [flux_io.max_snapshot_filename_length]u8 = undefined;
    const pvd_name = std.fmt.bufPrint(&pvd_buf, "{s}.pvd", .{base_name}) catch
        return error.FilenameTooLong;

    var dir = try std_fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(pvd_name, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn ensureDir(path: []const u8) !void {
    std.fs.cwd().makeDir(path) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };
}

fn canonicalEdge(a: u32, b: u32) [2]u32 {
    return if (a < b) .{ a, b } else .{ b, a };
}

fn signedTriangleArea(a: [2]f64, b: [2]f64, c: [2]f64) f64 {
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
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

fn transpose2(matrix: [2][2]f64) [2][2]f64 {
    return .{
        .{ matrix[0][0], matrix[1][0] },
        .{ matrix[0][1], matrix[1][1] },
    };
}

fn mulMat2(left: [2][2]f64, right: [2][2]f64) [2][2]f64 {
    return .{
        .{
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        },
        .{
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        },
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

const projection_candidates = [_][3]f64{
    .{ 1.0, 2.0, 3.0 },
    .{ 1.0, 3.0, 2.0 },
    .{ 2.0, 1.0, 3.0 },
    .{ 2.0, 3.0, 1.0 },
    .{ 3.0, 1.0, 2.0 },
    .{ 3.0, 2.0, 1.0 },
    .{ -1.0, 2.0, 3.0 },
    .{ 1.0, -2.0, 3.0 },
    .{ 1.0, 2.0, -3.0 },
    .{ -2.0, 1.0, 3.0 },
    .{ 2.0, -1.0, 3.0 },
    .{ 2.0, 1.0, -3.0 },
    .{ -3.0, 1.0, 2.0 },
    .{ 3.0, -1.0, 2.0 },
    .{ 3.0, 1.0, -2.0 },
    .{ 1.0, 1.0, 1.0 },
    .{ 1.0, 1.0, -1.0 },
    .{ 1.0, -1.0, 1.0 },
    .{ -1.0, 1.0, 1.0 },
};

test "surface diffusion error decreases under sphere refinement" {
    const allocator = testing.allocator;
    const refinements = [_]u32{ 1, 2, 3 };

    const results = try runConvergenceStudy(allocator, &refinements);
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

    const result = try run(allocator, .{
        .refinement = 0,
        .steps = 2,
        .frames = 0,
        .output_dir = "output/diffusion-surface-test",
    }, log_buffer.writer(allocator));

    try testing.expectEqual(@as(u32, 2), result.steps);
    try testing.expect(result.elapsed_s >= 0.0);
}
