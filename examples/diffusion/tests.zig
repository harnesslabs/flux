const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const diffusion = @import("root.zig");

test "surface diffusion error decreases under sphere refinement" {
    const allocator = testing.allocator;
    const refinements = [_]u32{ 1, 2, 3 };

    const results = try diffusion.runConvergenceStudy(.sphere, allocator, &refinements);
    defer allocator.free(results);

    try testing.expectEqual(refinements.len, results.len);
    try testing.expect(results[1].l2_error < results[0].l2_error);
    try testing.expect(results[2].l2_error < results[1].l2_error);
    try testing.expect(results[2].l2_error < 1e-3);
}

test "heat convergence study reaches second-order spatial rate" {
    const allocator = testing.allocator;
    const grids = [_]u32{ 8, 16, 32 };
    const results = try diffusion.runConvergenceStudy(.plane, allocator, &grids);
    defer allocator.free(results);
    var errors = [_]f64{0} ** grids.len;
    for (results, 0..) |result, idx| {
        errors[idx] = result.l2_error;
    }
    const rates = try flux.integrators.evolution.empiricalRates(allocator, &errors, 2.0);
    defer allocator.free(rates);

    try testing.expectEqual(grids.len, results.len);
    try testing.expectEqual(results.len - 1, rates.len);
    for (rates) |rate| {
        try testing.expect(rate > 1.75);
    }
}

test {
    testing.refAllDeclsRecursive(diffusion);
}
