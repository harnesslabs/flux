const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

pub const Mesh2D = flux.Mesh(2, 2);
pub const VertexField = flux.Cochain(Mesh2D, 0, flux.Primal);

pub const Config = struct {
    grid: u32 = 8,
    steps: u32 = 8,
    domain: f64 = 1.0,
    dt_scale: f64 = 0.1,
    output_dir: []const u8 = "output/heat",
    frames: u32 = 4,

    pub fn dt(self: Config) f64 {
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

pub fn run(
    allocator: std.mem.Allocator,
    config: Config,
    writer: anytype,
) !RunResult {
    _ = allocator;
    _ = writer;
    _ = config;
    return error.NotYetImplemented;
}

pub fn runConvergenceStudy(
    allocator: std.mem.Allocator,
    grids: []const u32,
) ![]ConvergenceResult {
    _ = allocator;
    _ = grids;
    return error.NotYetImplemented;
}

test "heat zero initial data stays zero under backward Euler" {
    const allocator = testing.allocator;
    const grids = [_]u32{4};
    const results = try runConvergenceStudy(allocator, &grids);
    defer allocator.free(results);

    try testing.expectEqual(@as(usize, 1), results.len);
    try testing.expectApproxEqAbs(@as(f64, 0.0), results[0].l2_error, 1e-12);
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
