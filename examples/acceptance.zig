//! M3 milestone capstone — end-to-end conservation invariants in one place.

const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

const maxwell = @import("maxwell");
const euler = @import("euler");
const diffusion = @import("diffusion");

const acceptance_steps: u32 = 10;

test "M3 acceptance: Maxwell 3D enforces ∇·B = 0 to machine precision over a short cavity run" {
    const allocator = testing.allocator;
    const Mesh = flux.topology.Mesh(3, 3);
    const Maxwell = maxwell.system_api.Maxwell(3, Mesh);
    const dt: f64 = 0.0025;

    var mesh = try Mesh.cartesian(allocator, .{ 2, 2, 2 }, .{ 1.0, 1.0, 1.0 });
    defer mesh.deinit(allocator);

    var system = try Maxwell.cavity(allocator, &mesh, .{
        .extents = .{ 1.0, 1.0, 1.0 },
        .time_step = dt,
        .boundary = .pec,
    });
    defer system.deinit(allocator);

    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try flux.integrators.Leapfrog(*Maxwell).advance(allocator, &system, dt);

        var derivative = try (try system.state.dec_operators.exteriorDerivative(flux.forms.Primal, 2)).apply(allocator, system.state.B);
        defer derivative.deinit(allocator);

        var divergence_sq: f64 = 0.0;
        for (derivative.values) |value| {
            divergence_sq += value * value;
        }
        try testing.expectApproxEqAbs(@as(f64, 0.0), @sqrt(divergence_sq), 1e-12);
    }
}

test "M3 acceptance: Euler 2D conserves total circulation to machine precision over a short dipole run" {
    const allocator = testing.allocator;
    const Mesh = flux.topology.Mesh(2, 2);
    const Euler = euler.system_api.Euler(2, Mesh);
    const dt: f64 = 0.1 * (1.0 / 16.0);

    var mesh = try Mesh.cartesian(allocator, .{ 16, 16 }, .{ 1.0, 1.0 });
    defer mesh.deinit(allocator);

    var system = try Euler.dipole(allocator, &mesh);
    defer system.deinit(allocator);

    const circulation_initial = try system.measurement(allocator, .circulation);
    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try Euler.Explicit.advance(allocator, &system, dt);
    }

    const circulation_final = try system.measurement(allocator, .circulation);
    try testing.expectApproxEqAbs(circulation_initial, circulation_final, 1e-12);
}

test "M3 acceptance: Euler 3D conserves helicity to machine precision over a short reference-mode run" {
    const allocator = testing.allocator;
    const Mesh = flux.topology.Mesh(3, 3);
    const Euler = euler.system_api.Euler(3, Mesh);

    var mesh = try Mesh.cartesian(allocator, .{ 2, 2, 2 }, .{ 1.0, 1.0, 1.0 });
    defer mesh.deinit(allocator);

    var system = try Euler.reference(allocator, &mesh);
    defer system.deinit(allocator);

    const helicity_initial = try system.measurement(allocator, .helicity);
    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try Euler.Explicit.advance(allocator, &system, 0.01);
    }

    const helicity_final = try system.measurement(allocator, .helicity);
    try testing.expectApproxEqAbs(helicity_initial, helicity_final, 1e-12);
}

test "M3 acceptance: heat equation reaches expected spatial convergence rate" {
    const allocator = testing.allocator;
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const plane8 = try diffusion.runPlane(allocator, .{
        .counts = .{ 8, 8 },
        .steps = 8,
        .extents = .{ 1.0, 1.0 },
        .snapshot_cadence = .disabled,
    }, stderr);
    const plane16 = try diffusion.runPlane(allocator, .{
        .counts = .{ 16, 16 },
        .steps = 32,
        .extents = .{ 1.0, 1.0 },
        .snapshot_cadence = .disabled,
    }, stderr);

    const errors = [_]f64{ plane8.summary.l2_error_final, plane16.summary.l2_error_final };
    const rates = try flux.evolution.reference.empiricalRates(allocator, &errors, 2.0);
    defer allocator.free(rates);

    try testing.expectEqual(@as(usize, 1), rates.len);
    try testing.expect(rates[0] > 1.75);
}

test "M3 acceptance: diffusion-surface solution matches the analytic eigenmode under refinement" {
    const allocator = testing.allocator;
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const coarse = try diffusion.runSphere(allocator, .{
        .refinement = 1,
        .steps = 8,
        .snapshot_cadence = .disabled,
    }, stderr);
    const fine = try diffusion.runSphere(allocator, .{
        .refinement = 2,
        .steps = 8,
        .snapshot_cadence = .disabled,
    }, stderr);

    try testing.expect(fine.summary.l2_error_final < coarse.summary.l2_error_final);
    try testing.expect(fine.summary.l2_error_final < 5e-3);
}
