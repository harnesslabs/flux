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

pub fn faradayStep(_: std.mem.Allocator, _: anytype, _: f64) !void {
    return error.NotImplemented;
}

pub fn ampereStep(_: std.mem.Allocator, _: anytype, _: f64) !void {
    return error.NotImplemented;
}

pub fn applyPecBoundary(_: std.mem.Allocator, _: anytype) !void {
    return error.NotImplemented;
}

pub fn leapfrogStep(_: std.mem.Allocator, _: anytype, _: f64) !void {
    return error.NotImplemented;
}

pub fn runSimulation(_: std.mem.Allocator, _: *MaxwellState3D, _: Config) !void {
    return error.NotImplemented;
}

pub fn writeSnapshot(_: std.mem.Allocator, _: anytype, _: *const MaxwellState3D) !void {
    return error.NotImplemented;
}

pub fn runCli() !void {
    return error.NotImplemented;
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
