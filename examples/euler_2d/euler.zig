const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const std_fs = std.fs;

pub const Mesh2D = flux.Mesh(2, 2);
pub const VertexVorticity = flux.Cochain(Mesh2D, 0, flux.Primal);
pub const FaceVorticity = flux.Cochain(Mesh2D, 2, flux.Primal);
pub const FaceTracer = flux.Cochain(Mesh2D, 2, flux.Primal);
const poisson = flux.operators.poisson;
const flux_io = flux.io;

const Vec2 = [2]f64;

pub const Demo = enum {
    gaussian,
    dipole,
};

pub const Config = struct {
    demo: Demo = .gaussian,
    grid: u32 = 16,
    steps: u32 = 1000,
    domain: f64 = 1.0,
    cfl: f64 = 0.1,
    output_dir: []const u8 = "output/euler_2d",
    frames: u32 = 50,

    pub fn dt(self: Config) f64 {
        return stableDt(self);
    }

    pub fn outputInterval(self: Config) u32 {
        if (self.frames == 0) return 0;
        return @max(1, self.steps / self.frames);
    }
};

pub const RunResult = struct {
    elapsed_s: f64,
    circulation_initial: f64,
    circulation_final: f64,
    snapshot_count: u32,
};

pub const State = struct {
    const EdgeAdjacency = struct {
        faces: [2]?u32 = .{ null, null },
    };

    mesh: *const Mesh2D,
    operators: *flux.OperatorContext(Mesh2D),
    stream_function: VertexVorticity,
    vorticity: FaceVorticity,
    tracer: FaceTracer,
    face_velocity: []Vec2,
    edge_adjacency: []EdgeAdjacency,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh2D) !State {
        var stream_function = try VertexVorticity.init(allocator, mesh);
        errdefer stream_function.deinit(allocator);

        var vorticity = try FaceVorticity.init(allocator, mesh);
        errdefer vorticity.deinit(allocator);

        var tracer = try FaceTracer.init(allocator, mesh);
        errdefer tracer.deinit(allocator);

        const operators = try flux.OperatorContext(Mesh2D).init(allocator, mesh);
        errdefer operators.deinit();
        try operators.withLaplacian(0);

        const face_velocity = try allocator.alloc(Vec2, mesh.num_faces());
        errdefer allocator.free(face_velocity);
        @memset(face_velocity, .{ 0.0, 0.0 });

        const edge_adjacency = try buildEdgeAdjacency(allocator, mesh);
        errdefer allocator.free(edge_adjacency);

        return .{
            .mesh = mesh,
            .operators = operators,
            .stream_function = stream_function,
            .vorticity = vorticity,
            .tracer = tracer,
            .face_velocity = face_velocity,
            .edge_adjacency = edge_adjacency,
        };
    }

    pub fn deinit(self: *State, allocator: std.mem.Allocator) void {
        allocator.free(self.edge_adjacency);
        allocator.free(self.face_velocity);
        self.tracer.deinit(allocator);
        self.vorticity.deinit(allocator);
        self.stream_function.deinit(allocator);
        self.operators.deinit();
    }
};

pub fn initializeGaussianVortex(state: *State) void {
    const face_centers = state.mesh.simplices(2).items(.barycenter);
    const sigma = 0.12;

    for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
        omega.* = gaussianBlob(center, .{ 0.5, 0.5 }, sigma);
        tracer.* = tracerStripe(center);
    }
}

pub fn initializeVortexDipole(state: *State) void {
    const face_centers = state.mesh.simplices(2).items(.barycenter);
    const sigma = 0.075;

    for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
        const positive = gaussianBlob(center, .{ 0.38, 0.50 }, sigma);
        const negative = gaussianBlob(center, .{ 0.62, 0.50 }, sigma);
        omega.* = 1.4 * (positive - negative);
        tracer.* = tracerStripe(center);
    }
}

pub fn totalCirculation(state: *const State) f64 {
    const areas = state.mesh.simplices(2).items(.volume);
    var circulation: f64 = 0.0;
    for (state.vorticity.values, areas) |omega, area| {
        circulation += omega * area;
    }
    return circulation;
}

pub fn step(
    allocator: std.mem.Allocator,
    state: *State,
    dt: f64,
) !void {
    try refreshDerivedFields(allocator, state);
    try advectVorticity(allocator, state, dt);
    try advectScalar(allocator, state, state.tracer.values, dt);
    try refreshDerivedFields(allocator, state);
}

fn stableDt(config: Config) f64 {
    return config.cfl * (config.domain / @as(f64, @floatFromInt(config.grid)));
}

pub fn run(
    allocator: std.mem.Allocator,
    config: Config,
    writer: anytype,
) !RunResult {
    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    switch (config.demo) {
        .gaussian => initializeGaussianVortex(&state),
        .dipole => initializeVortexDipole(&state),
    }
    try refreshDerivedFields(allocator, &state);

    const circulation_initial = totalCirculation(&state);
    const dt = config.dt();
    const interval = config.outputInterval();
    const has_output = interval > 0;
    const snapshot_capacity: u32 = if (has_output) (config.steps / interval) +| 1 else 0;

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

    const start_ns = std.time.nanoTimestamp();
    var snapshot_count: u32 = 0;
    for (0..config.steps) |step_idx| {
        try step(allocator, &state, dt);

        if (has_output and (step_idx + 1) % interval == 0) {
            const filename = flux_io.snapshot_filename(
                &filename_bufs[snapshot_count],
                baseName(config.demo),
                snapshot_count,
            );
            try writeSnapshot(allocator, config.output_dir, filename, &state);
            pvd_entries[snapshot_count] = .{
                .timestep = @as(f64, @floatFromInt(step_idx + 1)) * dt,
                .filename = filename,
            };
            snapshot_count += 1;
        }
    }
    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

    if (has_output and snapshot_count > 0) {
        try writePvd(allocator, config.output_dir, baseName(config.demo), pvd_entries[0..snapshot_count]);
    }

    const circulation_final = totalCirculation(&state);
    try writer.print(
        "euler_2d[{s}]: grid={d} steps={d} dt={d:.6} circulation={d:.12} -> {d:.12}\n",
        .{ @tagName(config.demo), config.grid, config.steps, dt, circulation_initial, circulation_final },
    );

    return .{
        .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
        .circulation_initial = circulation_initial,
        .circulation_final = circulation_final,
        .snapshot_count = snapshot_count,
    };
}

fn baseName(demo: Demo) []const u8 {
    return switch (demo) {
        .gaussian => "euler_2d",
        .dipole => "euler_dipole",
    };
}

fn buildEdgeAdjacency(allocator: std.mem.Allocator, mesh: *const Mesh2D) ![]State.EdgeAdjacency {
    const adjacency = try allocator.alloc(State.EdgeAdjacency, mesh.num_edges());
    @memset(adjacency, .{ .faces = .{ null, null } });

    const boundary_2 = mesh.boundary(2);
    for (0..mesh.num_faces()) |face_idx_usize| {
        const face_idx: u32 = @intCast(face_idx_usize);
        const row = boundary_2.row(face_idx);
        for (row.cols) |edge_idx| {
            const edge = &adjacency[edge_idx];
            if (edge.faces[0] == null) {
                edge.faces[0] = face_idx;
                continue;
            }
            std.debug.assert(edge.faces[1] == null);
            edge.faces[1] = face_idx;
        }
    }

    return adjacency;
}

fn refreshDerivedFields(allocator: std.mem.Allocator, state: *State) !void {
    try recoverStreamFunction(allocator, state);
    reconstructFaceVelocity(state);
}

fn recoverStreamFunction(allocator: std.mem.Allocator, state: *State) !void {
    const vertex_count = state.mesh.num_vertices();
    const forcing = try allocator.alloc(f64, vertex_count);
    defer allocator.free(forcing);
    const weights = try allocator.alloc(f64, vertex_count);
    defer allocator.free(weights);
    const boundary = try allocator.alloc(f64, vertex_count);
    defer allocator.free(boundary);

    @memset(forcing, 0.0);
    @memset(weights, 0.0);
    @memset(boundary, 0.0);

    const face_vertices = state.mesh.simplices(2).items(.vertices);
    const face_areas = state.mesh.simplices(2).items(.volume);

    for (face_vertices, face_areas, state.vorticity.values) |verts, area, omega| {
        const lumped_area = area / 3.0;
        for (verts) |vertex_idx| {
            forcing[vertex_idx] += omega * lumped_area;
            weights[vertex_idx] += lumped_area;
        }
    }

    for (forcing, weights) |*forcing_value, weight| {
        if (weight == 0.0) continue;
        forcing_value.* /= weight;
    }

    var solve = try poisson.solve_zero_form_dirichlet(Mesh2D, allocator, state.operators, forcing, boundary, .{});
    defer solve.deinit(allocator);
    @memcpy(state.stream_function.values, solve.solution);
}

fn reconstructFaceVelocity(state: *State) void {
    const coords = state.mesh.vertices.slice().items(.coords);
    const face_vertices = state.mesh.simplices(2).items(.vertices);

    for (state.face_velocity, face_vertices) |*velocity, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const p2 = coords[verts[2]];

        const psi0 = state.stream_function.values[verts[0]];
        const psi1 = state.stream_function.values[verts[1]];
        const psi2 = state.stream_function.values[verts[2]];

        const det = (p1[0] - p0[0]) * (p2[1] - p0[1]) -
            (p2[0] - p0[0]) * (p1[1] - p0[1]);
        std.debug.assert(det > 0.0);

        const grad_x = (psi0 * (p1[1] - p2[1]) +
            psi1 * (p2[1] - p0[1]) +
            psi2 * (p0[1] - p1[1])) / det;
        const grad_y = (psi0 * (p2[0] - p1[0]) +
            psi1 * (p0[0] - p2[0]) +
            psi2 * (p1[0] - p0[0])) / det;

        velocity.* = .{ grad_y, -grad_x };
    }
}

fn advectVorticity(
    allocator: std.mem.Allocator,
    state: *State,
    dt: f64,
) !void {
    try advectScalar(allocator, state, state.vorticity.values, dt);
}

fn advectScalar(
    allocator: std.mem.Allocator,
    state: *const State,
    values: []f64,
    dt: f64,
) !void {
    const face_count = state.mesh.num_faces();
    const face_areas = state.mesh.simplices(2).items(.volume);
    const face_centers = state.mesh.simplices(2).items(.barycenter);
    const edge_vertices = state.mesh.simplices(1).items(.vertices);
    const edge_lengths = state.mesh.simplices(1).items(.volume);
    const coords = state.mesh.vertices.slice().items(.coords);

    const mass_old = try allocator.alloc(f64, face_count);
    defer allocator.free(mass_old);
    const mass_new = try allocator.alloc(f64, face_count);
    defer allocator.free(mass_new);

    for (mass_old, mass_new, values, face_areas) |*old, *new, value, area| {
        old.* = value * area;
        new.* = old.*;
    }

    for (state.edge_adjacency, 0..) |adjacency, edge_idx_usize| {
        const face_left = adjacency.faces[0] orelse continue;
        const face_right = adjacency.faces[1] orelse continue;

        const edge_idx: u32 = @intCast(edge_idx_usize);
        const edge = edge_vertices[edge_idx];
        const p0 = coords[edge[0]];
        const p1 = coords[edge[1]];
        const tangent = .{ p1[0] - p0[0], p1[1] - p0[1] };

        var normal = .{ tangent[1] / edge_lengths[edge_idx], -tangent[0] / edge_lengths[edge_idx] };
        const center_delta = .{
            face_centers[face_right][0] - face_centers[face_left][0],
            face_centers[face_right][1] - face_centers[face_left][1],
        };
        if (dot(normal, center_delta) < 0.0) {
            normal = .{ -normal[0], -normal[1] };
        }

        const velocity_left = state.face_velocity[face_left];
        const velocity_right = state.face_velocity[face_right];
        const edge_velocity = .{
            0.5 * (velocity_left[0] + velocity_right[0]),
            0.5 * (velocity_left[1] + velocity_right[1]),
        };

        const edge_flux = dt * edge_lengths[edge_idx] * dot(edge_velocity, normal);
        const transfer = if (edge_flux >= 0.0)
            edge_flux * values[face_left]
        else
            edge_flux * values[face_right];

        mass_new[face_left] -= transfer;
        mass_new[face_right] += transfer;
    }

    for (values, mass_new, face_areas) |*value, mass, area| {
        value.* = mass / area;
    }
}

fn dot(lhs: Vec2, rhs: Vec2) f64 {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1];
}

fn writeSnapshot(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    filename: []const u8,
    state: *const State,
) !void {
    var velocity_values = try allocator.alloc(f64, state.mesh.num_faces() * 3);
    defer allocator.free(velocity_values);

    for (state.face_velocity, 0..) |velocity, face_idx| {
        const base = 3 * face_idx;
        velocity_values[base + 0] = velocity[0];
        velocity_values[base + 1] = velocity[1];
        velocity_values[base + 2] = 0.0;
    }

    const point_data = [_]flux_io.DataArraySlice{
        .{ .name = "stream_function", .values = state.stream_function.values },
    };
    const cell_data = [_]flux_io.DataArraySlice{
        .{ .name = "vorticity", .values = state.vorticity.values },
        .{ .name = "tracer", .values = state.tracer.values },
        .{ .name = "velocity", .values = velocity_values, .num_components = 3 },
    };

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);
    try flux_io.write(output.writer(allocator), 2, 2, state.mesh.*, &point_data, &cell_data);

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

fn gaussianBlob(point: Vec2, center: Vec2, sigma: f64) f64 {
    const dx = point[0] - center[0];
    const dy = point[1] - center[1];
    return std.math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

fn tracerStripe(point: Vec2) f64 {
    const y = point[1] - 0.5;
    const stripe = std.math.exp(-(y * y) / (2.0 * 0.035 * 0.035));
    const x_modulation = 0.5 * (1.0 + std.math.cos(10.0 * std.math.pi * point[0]));
    return stripe * x_modulation;
}

test "Euler 2D state allocates stream function on vertices and vorticity on faces" {
    const allocator = testing.allocator;

    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_vertices()), state.stream_function.values.len);
    try testing.expectEqual(@as(usize, mesh.num_faces()), state.vorticity.values.len);
    try testing.expectEqual(@as(usize, mesh.num_faces()), state.tracer.values.len);
}

test "Euler 2D dipole initialization contains both circulation signs" {
    const allocator = testing.allocator;

    var mesh = try Mesh2D.uniform_grid(allocator, 10, 10, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeVortexDipole(&state);

    var min_value: f64 = std.math.inf(f64);
    var max_value: f64 = -std.math.inf(f64);
    for (state.vorticity.values) |omega| {
        min_value = @min(min_value, omega);
        max_value = @max(max_value, omega);
    }

    try testing.expect(min_value < 0.0);
    try testing.expect(max_value > 0.0);
}

test "Euler 2D circulation is conserved over 1000 explicit steps" {
    const allocator = testing.allocator;
    const config = Config{
        .grid = 8,
        .steps = 1000,
        .domain = 1.0,
        .cfl = 0.02,
    };

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    initializeGaussianVortex(&state);

    const circulation_initial = totalCirculation(&state);
    const dt = stableDt(config);

    for (0..config.steps) |_| {
        try step(allocator, &state, dt);
    }

    const circulation_final = totalCirculation(&state);
    try testing.expectApproxEqAbs(circulation_initial, circulation_final, 1e-12);
}
