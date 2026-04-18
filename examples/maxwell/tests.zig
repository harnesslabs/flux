const std = @import("std");
const testing = std.testing;
const root = @import("root.zig");
const reference = @import("reference.zig");

test "TE₁₀ cavity: eigenvalue error reduces ≥3× when grid halves (O(h²) convergence)" {
    const allocator = testing.allocator;
    const domain_length: f64 = 1.0;
    const analytical_omega_sq = std.math.pi * std.math.pi / (domain_length * domain_length);

    const omega_sq_coarse = try reference.compute_te10_eigenvalue(allocator, 8, domain_length);
    const omega_sq_fine = try reference.compute_te10_eigenvalue(allocator, 16, domain_length);

    const error_coarse = @abs(omega_sq_coarse - analytical_omega_sq) / analytical_omega_sq;
    const error_fine = @abs(omega_sq_fine - analytical_omega_sq) / analytical_omega_sq;

    try testing.expect(omega_sq_coarse > 0.0);
    try testing.expect(omega_sq_fine > 0.0);
    try testing.expect(error_fine < error_coarse);
    try testing.expect(error_coarse / error_fine >= 3.0);
}

test "convergence: energy drift bounded for TE₁₀ cavity" {
    const allocator = testing.allocator;

    const final_time: f64 = 0.5;
    const drift_coarse = try reference.run_cavity_energy_drift(allocator, 8, final_time);
    const drift_fine = try reference.run_cavity_energy_drift(allocator, 16, final_time);

    try testing.expect(drift_coarse < 0.02);
    try testing.expect(drift_fine < 0.02);
}

test "TM₁₁₀ cavity: eigenvalue error decreases under 3D refinement" {
    const allocator = testing.allocator;
    const analytical = 2.0 * std.math.pi * std.math.pi;

    const coarse = try reference.compute_tm110_eigenvalue(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    const fine = try reference.compute_tm110_eigenvalue(allocator, 4, 4, 4, 1.0, 1.0, 1.0);

    const error_coarse = @abs(coarse - analytical) / analytical;
    const error_fine = @abs(fine - analytical) / analytical;

    try testing.expect(error_fine < error_coarse);
    try testing.expect(error_coarse / error_fine >= 1.5);
}

test "TE₁₀ cavity reference measurements shrink under 2D refinement" {
    const allocator = testing.allocator;
    var writer_buffer: [1024]u8 = undefined;
    var writer_impl = std.Io.Writer.fixed(&writer_buffer);
    const writer = &writer_impl.interface;

    const coarse = try root.runCavityWithReference(2, allocator, .{
        .grid = 8,
        .steps = 12,
        .frames = 0,
        .output_dir = "output/maxwell-test-reference-coarse",
    }, writer);
    const fine = try root.runCavityWithReference(2, allocator, .{
        .grid = 16,
        .steps = 24,
        .frames = 0,
        .output_dir = "output/maxwell-test-reference-fine",
    }, writer);

    try testing.expect(fine.electric_l2_final < coarse.electric_l2_final);
    try testing.expect(fine.magnetic_l2_final < coarse.magnetic_l2_final);
}

test {
    testing.refAllDeclsRecursive(root);
}
