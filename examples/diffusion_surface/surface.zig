const std = @import("std");
const testing = std.testing;

pub const Config = struct {
    refinement: u32 = 0,
    steps: u32 = 8,
    dt_scale: f64 = 0.1,
    final_time: f64 = 0.05,
    output_dir: []const u8 = "output/diffusion_surface",
    frames: u32 = 4,
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

pub fn run(
    allocator: std.mem.Allocator,
    config: Config,
    writer: anytype,
) !RunResult {
    _ = allocator;
    _ = config;
    _ = writer;
    return error.NotYetImplemented;
}

pub fn runConvergenceStudy(
    allocator: std.mem.Allocator,
    refinements: []const u32,
) ![]ConvergenceResult {
    _ = allocator;
    _ = refinements;
    return error.NotYetImplemented;
}

test "surface diffusion error decreases under sphere refinement" {
    const allocator = testing.allocator;
    const refinements = [_]u32{ 0, 1, 2 };

    const results = try runConvergenceStudy(allocator, &refinements);
    defer allocator.free(results);

    try testing.expectEqual(refinements.len, results.len);
    try testing.expect(results[1].l2_error < results[0].l2_error);
    try testing.expect(results[2].l2_error < results[1].l2_error);
    try testing.expect(results[2].l2_error < 8e-2);
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
