const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

const cochain = flux.forms;
const topology = flux.topology;

pub const Mesh3D = topology.Mesh(3, 3);

pub const Config = struct {
    steps: u32 = 1000,
    nx: u32 = 2,
    ny: u32 = 2,
    nz: u32 = 2,
    width: f64 = 1.0,
    height: f64 = 1.0,
    depth: f64 = 1.0,
    dt: f64 = 0.01,
    output_dir: ?[]const u8 = null,
    output_interval: u32 = 0,
};

pub fn State(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const OneForm = cochain.Cochain(MeshType, 1, cochain.Primal);
        pub const TwoForm = cochain.Cochain(MeshType, 2, cochain.Primal);

        E: OneForm,
        B: TwoForm,
        J: OneForm,
        mesh: *const MeshType,
        operators: *flux.OperatorContext(MeshType),
        timestep: u64,

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var electric = try OneForm.init(allocator, mesh);
            errdefer electric.deinit(allocator);

            var magnetic = try TwoForm.init(allocator, mesh);
            errdefer magnetic.deinit(allocator);

            var current = try OneForm.init(allocator, mesh);
            errdefer current.deinit(allocator);

            const operators = try flux.OperatorContext(MeshType).init(allocator, mesh);
            errdefer operators.deinit();

            try operators.withExteriorDerivative(cochain.Primal, 1);
            try operators.withExteriorDerivative(cochain.Primal, 2);
            try operators.withExteriorDerivative(cochain.Dual, 1);
            try operators.withHodgeStar(2);
            try operators.withHodgeStarInverse(1);

            return .{
                .E = electric,
                .B = magnetic,
                .J = current,
                .mesh = mesh,
                .operators = operators,
                .timestep = 0,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.operators.deinit();
            self.J.deinit(allocator);
            self.B.deinit(allocator);
            self.E.deinit(allocator);
        }
    };
}

pub const MaxwellState3D = State(Mesh3D);

pub fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    var derivative = try state.operators.exteriorDerivative(cochain.Primal, 1).apply(allocator, state.E);
    defer derivative.deinit(allocator);

    derivative.scale(dt);
    state.B.sub(derivative);
}

pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    var star_b = try state.operators.hodgeStar(2).apply(allocator, state.B);
    defer star_b.deinit(allocator);

    var derivative = try state.operators.exteriorDerivative(cochain.Dual, 1).apply(allocator, star_b);
    defer derivative.deinit(allocator);

    var curl_b = try state.operators.hodgeStarInverse(1).apply(allocator, derivative);
    defer curl_b.deinit(allocator);

    for (state.E.values, curl_b.values, state.J.values) |*electric, curl_value, current| {
        electric.* += dt * (curl_value - current);
    }
}

pub fn applyPecBoundary(_: std.mem.Allocator, state: anytype) !void {
    for (state.mesh.boundary_edges) |edge_idx| {
        state.E.values[edge_idx] = 0.0;
    }
}

pub fn leapfrogStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try faradayStep(allocator, state, dt);
    try ampereStep(allocator, state, dt);
    try applyPecBoundary(allocator, state);
    state.timestep += 1;
}

pub fn runSimulation(allocator: std.mem.Allocator, state: *MaxwellState3D, config: Config) !void {
    var step_index: u32 = 0;
    while (step_index < config.steps) : (step_index += 1) {
        try leapfrogStep(allocator, state, config.dt);
    }
}

fn projectEdgesToTets(
    allocator: std.mem.Allocator,
    mesh: *const Mesh3D,
    edge_values: []const f64,
) ![]f64 {
    const tet_count = mesh.num_tets();
    const output = try allocator.alloc(f64, tet_count);
    errdefer allocator.free(output);
    const touched = try allocator.alloc(bool, mesh.num_edges());
    defer allocator.free(touched);

    for (0..tet_count) |tet_idx_usize| {
        @memset(touched, false);
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        var count: u32 = 0;

        for (tet_row.cols) |face_idx| {
            const face_row = mesh.boundary(2).row(face_idx);
            for (face_row.cols) |edge_idx| {
                if (touched[edge_idx]) continue;
                touched[edge_idx] = true;
                sum += @abs(edge_values[edge_idx]);
                count += 1;
            }
        }

        std.debug.assert(count > 0);
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(count));
    }

    return output;
}

fn projectFacesToTets(
    allocator: std.mem.Allocator,
    mesh: *const Mesh3D,
    face_values: []const f64,
) ![]f64 {
    const tet_count = mesh.num_tets();
    const output = try allocator.alloc(f64, tet_count);
    errdefer allocator.free(output);

    for (0..tet_count) |tet_idx_usize| {
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        for (tet_row.cols) |face_idx| {
            sum += @abs(face_values[face_idx]);
        }
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(tet_row.cols.len));
    }

    return output;
}

pub fn writeSnapshot(
    allocator: std.mem.Allocator,
    writer: anytype,
    state: *const MaxwellState3D,
) !void {
    const e_projected = try projectEdgesToTets(allocator, state.mesh, state.E.values);
    defer allocator.free(e_projected);

    const b_projected = try projectFacesToTets(allocator, state.mesh, state.B.values);
    defer allocator.free(b_projected);

    const cell_data = [_]flux.io.DataArraySlice{
        .{ .name = "E_intensity", .values = e_projected },
        .{ .name = "B_flux", .values = b_projected },
    };

    try flux.io.write(
        writer,
        Mesh3D.embedding_dimension,
        Mesh3D.topological_dimension,
        state.mesh.*,
        &.{},
        &cell_data,
    );
}

pub fn runCli() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var mesh = try makeCavityMesh(allocator, .{});
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedClosedMagneticField(allocator, &state);
    try runSimulation(allocator, &state, .{});

    var stdout_buffer = std.ArrayListUnmanaged(u8){};
    defer stdout_buffer.deinit(allocator);
    try writeSnapshot(allocator, stdout_buffer.writer(allocator), &state);
    const stdout = (std.fs.File{ .handle = std.posix.STDOUT_FILENO }).deprecatedWriter();
    try stdout.print(
        "maxwell_3d completed {d} steps on {d} tetrahedra\n",
        .{ state.timestep, mesh.num_tets() },
    );
}

fn makeCavityMesh(allocator: std.mem.Allocator, config: Config) !Mesh3D {
    return Mesh3D.uniform_tetrahedral_grid(
        allocator,
        config.nx,
        config.ny,
        config.nz,
        config.width,
        config.height,
        config.depth,
    );
}

fn seedClosedMagneticField(allocator: std.mem.Allocator, state: *MaxwellState3D) !void {
    var potential = try MaxwellState3D.OneForm.init(allocator, state.mesh);
    defer potential.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0x92_3D_B_00);
    for (potential.values) |*value| {
        value.* = rng.random().float(f64) * 0.2 - 0.1;
    }

    var exact_flux = try state.operators.exteriorDerivative(cochain.Primal, 1).apply(allocator, potential);
    defer exact_flux.deinit(allocator);

    @memcpy(state.B.values, exact_flux.values);
}

fn divergenceNorm(allocator: std.mem.Allocator, state: *const MaxwellState3D) !f64 {
    var divergence = try state.operators.exteriorDerivative(cochain.Primal, 2).apply(allocator, state.B);
    defer divergence.deinit(allocator);
    return std.math.sqrt(divergence.norm_squared());
}

test "MaxwellState3D initializes field spaces for tetrahedral meshes" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 1, .ny = 1, .nz = 1 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(@as(u64, 0), state.timestep);
    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.E.values.len)));
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(state.B.values.len)));
    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.J.values.len)));

    comptime {
        try testing.expectEqual(1, MaxwellState3D.OneForm.degree);
        try testing.expectEqual(2, MaxwellState3D.TwoForm.degree);
    }
}

test "3D Maxwell preserves d₂B = 0 over 1000 source-free cavity steps" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 2, .ny = 2, .nz = 2, .dt = 0.0025 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedClosedMagneticField(allocator, &state);

    for (0..1000) |_| {
        try leapfrogStep(allocator, &state, 0.0025);

        const norm = try divergenceNorm(allocator, &state);
        try testing.expectApproxEqAbs(@as(f64, 0.0), norm, 1e-12);
    }
}

test "3D Maxwell PEC zeroes boundary edge circulation after a step" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 2, .ny = 1, .nz = 1, .dt = 0.0025 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    for (state.E.values, 0..) |*value, idx| {
        value.* = @as(f64, @floatFromInt(idx + 1));
    }

    try applyPecBoundary(allocator, &state);

    for (mesh.boundary_edges) |edge_idx| {
        try testing.expectEqual(@as(f64, 0.0), state.E.values[edge_idx]);
    }
}

test "3D Maxwell snapshot exports projected E and B data on tetrahedral cells" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 1, .ny = 1, .nz = 1 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedClosedMagneticField(allocator, &state);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try writeSnapshot(allocator, output.writer(allocator), &state);

    const xml = output.items;
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"E_intensity\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"B_flux\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "NumberOfCells=\"6\"") != null);
}
