const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

const cochain = flux.forms;
const observers = flux.operators.observers;
const poisson = flux.operators.poisson;
const wedge_product = flux.operators.wedge_product;

pub const Mesh3D = flux.Mesh(3, 3);
pub const Velocity = flux.Cochain(Mesh3D, 1, flux.Primal);
pub const Vorticity = flux.Cochain(Mesh3D, 2, flux.Primal);
pub const Potential = flux.Cochain(Mesh3D, 2, flux.Primal);

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

pub const RunResult = struct {
    elapsed_s: f64,
    helicity_initial: f64,
    helicity_final: f64,
    snapshot_count: u32,
};

pub const State = struct {
    mesh: *const Mesh3D,
    operators: *flux.OperatorContext(Mesh3D),
    velocity: Velocity,
    vorticity: Vorticity,
    boundary_velocity: []f64,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh3D) !State {
        _ = allocator;
        _ = mesh;
        @panic("not yet implemented");
    }

    pub fn deinit(self: *State, allocator: std.mem.Allocator) void {
        _ = self;
        _ = allocator;
        @panic("not yet implemented");
    }
};

fn selectVelocity(state: *const State) *const Velocity {
    return &state.velocity;
}

pub fn seedReferenceMode(allocator: std.mem.Allocator, state: *State) !void {
    _ = allocator;
    _ = state;
    @panic("not yet implemented");
}

pub fn recoverVelocityFromVorticity(allocator: std.mem.Allocator, state: *State) !void {
    _ = allocator;
    _ = state;
    @panic("not yet implemented");
}

pub fn step(allocator: std.mem.Allocator, state: *State, dt: f64) !void {
    _ = allocator;
    _ = state;
    _ = dt;
    @panic("not yet implemented");
}

pub fn run(allocator: std.mem.Allocator, config: Config, writer: anytype) !RunResult {
    _ = allocator;
    _ = config;
    _ = writer;
    @panic("not yet implemented");
}

pub fn runCli() !void {
    @panic("not yet implemented");
}

test "velocity recovery reproduces the seeded co-closed 1-form on a tetrahedral mesh" {
    const allocator = testing.allocator;

    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedReferenceMode(allocator, &state);

    var expected = try Velocity.init(allocator, &mesh);
    defer expected.deinit(allocator);
    @memcpy(expected.values, state.velocity.values);

    try recoverVelocityFromVorticity(allocator, &state);

    for (state.velocity.values, expected.values) |actual, reference| {
        try testing.expectApproxEqAbs(reference, actual, 1e-10);
    }
}

test "helicity is conserved over 1000 steps for the seeded reference mode" {
    const allocator = testing.allocator;

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const result = try run(allocator, .{
        .steps = 1000,
        .nx = 2,
        .ny = 2,
        .nz = 2,
        .dt = 0.01,
        .output_dir = null,
        .output_interval = 0,
    }, output.writer(allocator));

    try testing.expectEqual(@as(u32, 0), result.snapshot_count);
    try testing.expectApproxEqAbs(result.helicity_initial, result.helicity_final, 1e-12);
}

test "seeded reference mode has nonzero helicity and matches the observer" {
    const allocator = testing.allocator;

    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    try seedReferenceMode(allocator, &state);

    const Helicity = observers.HelicityObserver(State, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    const observed = try observer.evaluate(allocator, &state, 0);

    var density = try wedge_product.wedge(allocator, state.velocity, state.vorticity);
    defer density.deinit(allocator);

    var manual: f64 = 0.0;
    for (density.values) |value| {
        manual += value;
    }

    try testing.expect(@abs(manual) > 1e-9);
    try testing.expectApproxEqAbs(manual, observed, 1e-12);
}
