const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

pub const Mesh2D = flux.Mesh(2, 2);
pub const VertexVorticity = flux.Cochain(Mesh2D, 0, flux.Primal);
pub const FaceVorticity = flux.Cochain(Mesh2D, 2, flux.Primal);

pub const Config = struct {
    grid: u32 = 16,
    steps: u32 = 1000,
    domain: f64 = 1.0,
    cfl: f64 = 0.1,
};

pub const State = struct {
    mesh: *const Mesh2D,
    operators: *flux.OperatorContext(Mesh2D),
    stream_function: VertexVorticity,
    vorticity: FaceVorticity,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh2D) !State {
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

pub fn initializeGaussianVortex(state: *State) void {
    _ = state;
    @panic("not yet implemented");
}

pub fn totalCirculation(state: *const State) f64 {
    _ = state;
    @panic("not yet implemented");
}

pub fn step(
    allocator: std.mem.Allocator,
    state: *State,
    dt: f64,
) !void {
    _ = allocator;
    _ = state;
    _ = dt;
    @panic("not yet implemented");
}

fn stableDt(config: Config) f64 {
    return config.cfl * (config.domain / @as(f64, @floatFromInt(config.grid)));
}

test "Euler 2D state allocates stream function on vertices and vorticity on faces" {
    const allocator = testing.allocator;

    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_vertices()), state.stream_function.values.len);
    try testing.expectEqual(@as(usize, mesh.num_faces()), state.vorticity.values.len);
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
