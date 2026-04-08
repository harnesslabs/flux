//! M3 milestone capstone — end-to-end conservation invariants in one place.
//!
//! Each existing example already verifies its own invariant in depth (1000-step
//! Maxwell ∇·B run, 1000-step helicity run, 8/16/32 heat refinement study, …).
//! This file is the *capstone*: short, focused tests that mirror the M3
//! milestone acceptance criterion verbatim and run on every CI push so a single
//! failing target proves M3 is broken without having to scan five files.
//!
//! Why these tests are short and complementary, not redundant. The per-example
//! tests prove each invariant under heavy load; the capstone proves all five
//! invariants hold together at the smallest grid that still demonstrates the
//! property, so CI surfaces a regression in the *milestone* contract — not just
//! in one example — within seconds.

const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

// Each physics example is imported as a named module (declared in build.zig)
// so its tests do not get pulled into this compilation unit. Only the five
// acceptance tests defined in this file run when this target executes.
const maxwell_3d = @import("maxwell_3d");
const euler_2d = @import("euler_2d");
const euler_3d = @import("euler_3d");
const heat = @import("heat");
const diffusion_surface = @import("diffusion_surface");

/// Step count for the short capstone runs. The per-example deep tests run for
/// 1000 steps; the capstone runs the same physics for an order of magnitude
/// fewer steps because the invariants asserted here are *structural* — they
/// hold to machine precision at every step or not at all.
const acceptance_steps: u32 = 10;

test "M3 acceptance: Maxwell 3D enforces ∇·B = 0 to machine precision over a short cavity run" {
    const allocator = testing.allocator;

    const config = maxwell_3d.Config{
        .nx = 2,
        .ny = 2,
        .nz = 2,
        .dt = 0.0025,
        .steps = acceptance_steps,
    };

    var mesh = try maxwell_3d.makeCavityMesh(allocator, config);
    defer mesh.deinit(allocator);

    var state = try maxwell_3d.MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try maxwell_3d.seedTm110Mode(allocator, &state, config.dt, config.width, config.height);

    // The invariant must hold at every step, not just at the end — a leapfrog
    // step that *temporarily* introduces ∇·B is already a failure.
    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try maxwell_3d.leapfrogStep(allocator, &state, config.dt);

        const divergence = try maxwell_3d.divergenceNorm(allocator, &state);
        try testing.expectApproxEqAbs(@as(f64, 0.0), divergence, 1e-12);
    }
}

test "M3 acceptance: Euler 2D conserves total circulation to machine precision over a short dipole run" {
    const allocator = testing.allocator;

    var mesh = try euler_2d.Mesh2D.uniform_grid(allocator, 16, 16, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try euler_2d.State.init(allocator, &mesh);
    defer state.deinit(allocator);

    euler_2d.initializeVortexDipole(&state);

    const dt: f64 = 0.1 * (1.0 / 16.0);
    const circulation_initial = euler_2d.totalCirculation(&state);

    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try euler_2d.step(allocator, &state, dt);
    }

    const circulation_final = euler_2d.totalCirculation(&state);
    try testing.expectApproxEqAbs(circulation_initial, circulation_final, 1e-12);
}

test "M3 acceptance: Euler 3D conserves helicity to machine precision over a short reference-mode run" {
    const allocator = testing.allocator;

    var mesh = try euler_3d.Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try euler_3d.State.init(allocator, &mesh);
    defer state.deinit(allocator);

    try euler_3d.seedReferenceMode(allocator, &state);
    const helicity_initial = try euler_3d.computeHelicity(allocator, &state);

    var step_idx: u32 = 0;
    while (step_idx < acceptance_steps) : (step_idx += 1) {
        try euler_3d.step(allocator, &state, 0.01);
    }

    const helicity_final = try euler_3d.computeHelicity(allocator, &state);
    try testing.expectApproxEqAbs(helicity_initial, helicity_final, 1e-12);
}

test "M3 acceptance: heat equation reaches expected spatial convergence rate" {
    const allocator = testing.allocator;

    // Two-grid pair is the smallest configuration that exhibits a measurable
    // rate. The deeper {8,16,32} sweep lives in the per-example test.
    const grids = [_]u32{ 8, 16 };
    const results = try heat.runConvergenceStudy(allocator, &grids);
    defer allocator.free(results);

    try testing.expectEqual(grids.len, results.len);
    const rate = std.math.log(f64, 2.0, results[0].l2_error / results[1].l2_error);
    // Backward-Euler heat solve is second order in space; require ≥ 1.75 to
    // accept asymptotic noise on the small-grid pair while still failing if
    // the operator stack regresses to first order.
    try testing.expect(rate > 1.75);
}

test "M3 acceptance: diffusion-surface solution matches the analytic eigenmode under refinement" {
    const allocator = testing.allocator;

    // Refinement levels 1 and 2 — the smallest pair that still shows error
    // strictly decreasing on the sphere; refinement 0 is too coarse for the
    // analytic eigenmode comparison.
    const refinements = [_]u32{ 1, 2 };
    const results = try diffusion_surface.runConvergenceStudy(allocator, &refinements);
    defer allocator.free(results);

    try testing.expectEqual(refinements.len, results.len);
    try testing.expect(results[1].l2_error < results[0].l2_error);
    // Loose absolute bound at refinement 2 — the per-example test asserts the
    // tighter < 1e-3 bound at refinement 3.
    try testing.expect(results[1].l2_error < 5e-3);
}
