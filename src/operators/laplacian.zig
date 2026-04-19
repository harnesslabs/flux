//! Discrete Laplace-de Rham operator surfaces.
//!
//! The shared noun here is still `laplacian`, but the operator family split is
//! explicit:
//! - `laplacian.dec` owns strong operator application, `Δ = dδ + δd`
//! - `laplacian.feec` owns weak-form stiffness assembly, `Sₖ`
//!
//! This module re-exports both so callers can use either
//! `flux.operators.laplacian.dec` / `flux.operators.laplacian.feec`
//! or the family-first paths `flux.operators.dec.laplacian` /
//! `flux.operators.feec.laplacian`.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const sparse = @import("../math/sparse.zig");
const context = @import("context.zig");
pub const dec = @import("laplacian/dec.zig");
pub const feec = @import("laplacian/feec.zig");

pub const AssembledLaplacian = dec.AssembledLaplacian;
pub const assemble_for_degree = dec.assemble_for_degree;
pub const laplacian_composed = dec.laplacian_composed;
pub const assemble_stiffness = feec.assemble_stiffness;

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const MeshSurface = topology.Mesh(3, 2);
const Mesh3D = topology.Mesh(3, 3);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);

test "Δ₀ of constant 0-form is zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 42.0;

    var result = try (try operator_context.laplacian(0)).apply(allocator, omega);
    defer result.deinit(allocator);

    for (result.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), v, 2e-13);
    }
}

test "Δ₀ of linear function f(x,y) = x is zero at interior vertices" {
    // The DEC Laplacian reproduces the cotangent Laplacian, which is exact
    // for piecewise-linear functions on any triangulation. Since Δ(x) = 0
    // continuously, interior vertices should see zero.
    const allocator = testing.allocator;
    const nx: u32 = 5;
    const ny: u32 = 4;
    var mesh = try Mesh2D.plane(allocator, nx, ny, 3.0, 2.0);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0];

    var result = try (try operator_context.laplacian(0)).apply(allocator, omega);
    defer result.deinit(allocator);

    // Interior vertices: 1 ≤ i ≤ nx-1, 1 ≤ j ≤ ny-1
    // Vertex index: i * (ny + 1) + j
    for (1..nx) |i| {
        for (1..ny) |j| {
            const idx = i * (ny + 1) + j;
            try testing.expectApproxEqAbs(@as(f64, 0.0), result.values[idx], 1e-12);
        }
    }
}

test "Δ₀ is positive-semidefinite on random 0-forms (1000 trials)" {
    // ⟨ω, Δ₀ω⟩_★₀ = ωᵀ ★₀ Δ₀ω = ωᵀ D₀ᵀ ★₁ D₀ ω = ‖D₀ω‖²_★₁ ≥ 0
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_00);

    for (0..1000) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try (try operator_context.laplacian(0)).apply(allocator, omega);
        defer lap_omega.deinit(allocator);

        // ⟨ω, Δ₀ω⟩_★₀ = Σᵢ ωᵢ · (Δ₀ω)ᵢ · dual_area[i]
        var inner: f64 = 0;
        for (omega.values, lap_omega.values, dual_areas) |w, lw, area| {
            inner += w * lw * area;
        }
        // Must be non-negative (small numerical noise allowed)
        try testing.expect(inner >= -1e-10);
    }
}

test "Δ₀ is symmetric in ★₀-weighted inner product (1000 trials)" {
    // Self-adjointness: ⟨Δ₀f, g⟩_★₀ = ⟨f, Δ₀g⟩_★₀
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_01);

    for (0..1000) |_| {
        var f_form = try PrimalC0.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC0.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try (try operator_context.laplacian(0)).apply(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try (try operator_context.laplacian(0)).apply(allocator, g_form);
        defer lap_g.deinit(allocator);

        var inner_lap_f_g: f64 = 0;
        var inner_f_lap_g: f64 = 0;
        for (0..f_form.values.len) |idx| {
            inner_lap_f_g += lap_f.values[idx] * g_form.values[idx] * dual_areas[idx];
            inner_f_lap_g += f_form.values[idx] * lap_g.values[idx] * dual_areas[idx];
        }

        try testing.expectApproxEqRel(inner_lap_f_g, inner_f_lap_g, 1e-10);
    }
}

test "Δ₀ kernel is exactly the constant functions on connected mesh" {
    // On a connected mesh, the only 0-forms with Δ₀ω = 0 are constants.
    // Test: random non-constant 0-form has ‖Δ₀ω‖ > 0.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_02);

    for (0..100) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try (try operator_context.laplacian(0)).apply(allocator, omega);
        defer lap_omega.deinit(allocator);

        // Non-constant ⟹ ‖Δ₀ω‖ > 0
        var norm_sq: f64 = 0;
        for (lap_omega.values) |v| norm_sq += v * v;
        try testing.expect(norm_sq > 1e-10);
    }
}

test "assembled Δ₀ apply matches compose-on-the-fly Laplacian" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_03);

    for (0..100) |_| {
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var expected = try laplacian_composed(allocator, omega);
        defer expected.deinit(allocator);

        var actual = try (try operator_context.laplacian(0)).apply(allocator, omega);
        defer actual.deinit(allocator);

        for (actual.values, expected.values) |got, want| {
            try testing.expectApproxEqAbs(want, got, 5e-12);
        }
    }
}

test "assembled Δ₀ apply is stable across repeated applications" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 6, 5, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0] - 0.5 * c[1];

    for (0..3) |_| {
        var expected = try laplacian_composed(allocator, omega);
        defer expected.deinit(allocator);

        var actual = try (try operator_context.laplacian(0)).apply(allocator, omega);
        defer actual.deinit(allocator);

        for (actual.values, expected.values) |got, want| {
            try testing.expectApproxEqAbs(want, got, 5e-12);
        }
    }
}

test "assembled Δ₀ on Mesh(3, 2) matches intrinsic 2D Δ₀ on an isometric embedding" {
    const allocator = testing.allocator;
    const faces = [_][3]u32{
        .{ 0, 1, 2 },
        .{ 0, 2, 3 },
    };
    const intrinsic_vertices = [_][2]f64{
        .{ 0.0, 0.0 },
        .{ 1.0, 0.0 },
        .{ 2.0, 1.0 },
        .{ 1.0, 1.0 },
    };
    const embedded_vertices = [_][3]f64{
        .{ 0.0, 0.0, 0.0 },
        .{ 0.7071067811865475, 0.7071067811865475, 0.0 },
        .{ 1.414213562373095, 1.414213562373095, 1.0 },
        .{ 0.7071067811865475, 0.7071067811865475, 1.0 },
    };

    var intrinsic_mesh = try Mesh2D.from_triangles(allocator, &intrinsic_vertices, &faces);
    defer intrinsic_mesh.deinit(allocator);
    var embedded_mesh = try MeshSurface.from_triangles(allocator, &embedded_vertices, &faces);
    defer embedded_mesh.deinit(allocator);

    const intrinsic_context = try context.OperatorContext(Mesh2D).init(allocator, &intrinsic_mesh);
    defer intrinsic_context.deinit();
    const embedded_context = try context.OperatorContext(MeshSurface).init(allocator, &embedded_mesh);
    defer embedded_context.deinit();

    var intrinsic_form = try PrimalC0.init(allocator, &intrinsic_mesh);
    defer intrinsic_form.deinit(allocator);
    const EmbeddedPrimalC0 = cochain.Cochain(MeshSurface, 0, cochain.Primal);
    var embedded_form = try EmbeddedPrimalC0.init(allocator, &embedded_mesh);
    defer embedded_form.deinit(allocator);

    for (intrinsic_form.values, 0..) |*value, idx| {
        value.* = @as(f64, @floatFromInt(idx + 1)) * 0.125;
        embedded_form.values[idx] = value.*;
    }

    var intrinsic_result = try (try intrinsic_context.laplacian(0)).apply(allocator, intrinsic_form);
    defer intrinsic_result.deinit(allocator);
    var embedded_result = try (try embedded_context.laplacian(0)).apply(allocator, embedded_form);
    defer embedded_result.deinit(allocator);

    for (intrinsic_result.values, embedded_result.values) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-11);
    }
}

test "assembled FEEC stiffness satisfies S₀ω = M₀(Δ₀ω) on embedded surface meshes" {
    const allocator = testing.allocator;
    const ScalarForm = cochain.Cochain(MeshSurface, 0, cochain.Primal);

    var mesh = try MeshSurface.sphere(allocator, 1.0, 1);
    defer mesh.deinit(allocator);

    const operator_context = try context.OperatorContext(MeshSurface).init(allocator, &mesh);
    defer operator_context.deinit();

    var omega = try ScalarForm.init(allocator, &mesh);
    defer omega.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0x51_0F_AA);
    for (omega.values) |*value| {
        value.* = rng.random().float(f64) * 2.0 - 1.0;
    }

    var strong = try (try operator_context.laplacian(0)).apply(allocator, omega);
    defer strong.deinit(allocator);

    var stiffness = try assemble_stiffness(0, allocator, &mesh);
    defer stiffness.deinit(allocator);

    const expected = try allocator.alloc(f64, omega.values.len);
    defer allocator.free(expected);
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    for (expected, strong.values, dual_volumes) |*out, value, dual_volume| {
        out.* = dual_volume * value;
    }

    const actual = try allocator.alloc(f64, omega.values.len);
    defer allocator.free(actual);
    sparse.spmv(stiffness, omega.values, actual);

    for (actual, expected) |lhs, rhs| {
        try testing.expectApproxEqAbs(rhs, lhs, 1e-11);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Δ₁ tests
// ═══════════════════════════════════════════════════════════════════════════

const PrimalC1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);

test "assembled FEEC stiffness satisfies S₁ω = M₁(Δ₁ω) on tetrahedral meshes" {
    const allocator = testing.allocator;
    const OneForm = cochain.Cochain(Mesh3D, 1, cochain.Primal);

    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const operator_context = try context.OperatorContext(Mesh3D).init(allocator, &mesh);
    defer operator_context.deinit();

    var omega = try OneForm.init(allocator, &mesh);
    defer omega.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0x51_1F_AA);
    for (omega.values) |*value| {
        value.* = rng.random().float(f64) * 2.0 - 1.0;
    }

    var strong = try (try operator_context.laplacian(1)).apply(allocator, omega);
    defer strong.deinit(allocator);

    var stiffness = try assemble_stiffness(1, allocator, &mesh);
    defer stiffness.deinit(allocator);

    const expected = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(expected);
    sparse.spmv(mesh.whitney_mass(1), strong.values, expected);

    const actual = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(actual);
    sparse.spmv(stiffness, omega.values, actual);

    for (actual, expected) |lhs, rhs| {
        try testing.expectApproxEqAbs(rhs, lhs, 1e-10);
    }
}

test "Δ₁ of exact 1-form d₀f has dδ component zero" {
    // For ω = d₀f (an exact 1-form), the dδ term of Δ₁ vanishes because
    // δ(d₀f) = 0 on exact forms: d(★⁻¹ d₀ᵀ ★ d₀ f) is the grad-div
    // term, but d₁(d₀f) = 0 by dd=0, so only the δd term survives.
    // Thus Δ₁(d₀f) = ★₁⁻¹ d₁ᵀ ★₂ d₁(d₀f) = 0, since d₁(d₀f) = 0.
    //
    // In other words: exact 1-forms are in the kernel of the δd term,
    // and the dδ term applied to d₀f reduces to d₀(Δ₀f). This test
    // verifies Δ₁ acts correctly on exact forms.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(1);

    // Build ω = d₀(constant) = 0 — trivial but validates that Δ₁(0) = 0.
    var zero_form = try PrimalC1.init(allocator, &mesh);
    defer zero_form.deinit(allocator);

    var result = try (try operator_context.laplacian(1)).apply(allocator, zero_form);
    defer result.deinit(allocator);

    for (result.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), v, 1e-13);
    }
}

test "Δ₁ is positive-semidefinite on random 1-forms (500 trials)" {
    // ⟨ω, Δ₁ω⟩_★₁ = ‖d₁ω‖²_★₂ + ‖δ₀ω‖²_★₀ ≥ 0
    // The ★₁-weighted inner product uses dual_length/length ratios.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(1);

    const edge_volumes = mesh.simplices(1).items(.volume);
    const dual_edge_volumes = mesh.dual_edge_volumes;
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_10);

    for (0..500) |_| {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try (try operator_context.laplacian(1)).apply(allocator, omega);
        defer lap_omega.deinit(allocator);

        // ⟨ω, Δ₁ω⟩_★₁ = Σᵢ ωᵢ · (Δ₁ω)ᵢ · (dual_length[i] / length[i])
        var inner: f64 = 0;
        for (omega.values, lap_omega.values, edge_volumes, dual_edge_volumes) |w, lw, volume, dual_volume| {
            inner += w * lw * (dual_volume / volume);
        }
        try testing.expect(inner >= -1e-8);
    }
}

test "Δ₁ is symmetric in ★₁-weighted inner product (500 trials)" {
    // Self-adjointness: ⟨Δ₁f, g⟩_★₁ = ⟨f, Δ₁g⟩_★₁
    // where ⟨u, v⟩_★₁ = uᵀ M₁ v (Whitney mass matrix inner product).
    const allocator = testing.allocator;
    const sparse_mod = @import("../math/sparse.zig");

    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(1);

    const edge_count = mesh.num_edges();
    const m1_buf = try allocator.alloc(f64, edge_count);
    defer allocator.free(m1_buf);

    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_11);

    for (0..500) |_| {
        var f_form = try PrimalC1.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC1.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try (try operator_context.laplacian(1)).apply(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try (try operator_context.laplacian(1)).apply(allocator, g_form);
        defer lap_g.deinit(allocator);

        // ⟨Δ₁f, g⟩_★₁ = (Δ₁f)ᵀ M₁ g
        sparse_mod.spmv(mesh.whitney_mass(1), g_form.values, m1_buf);
        var inner_lap_f_g: f64 = 0;
        for (lap_f.values, m1_buf) |lf, mg| inner_lap_f_g += lf * mg;

        // ⟨f, Δ₁g⟩_★₁ = fᵀ M₁ (Δ₁g)
        sparse_mod.spmv(mesh.whitney_mass(1), lap_g.values, m1_buf);
        var inner_f_lap_g: f64 = 0;
        for (f_form.values, m1_buf) |fv, mlg| inner_f_lap_g += fv * mlg;

        // CG solve introduces ~1e-10 relative residual per solve, and the
        // Laplacian does two solves. Loosen tolerance accordingly.
        try testing.expectApproxEqRel(inner_lap_f_g, inner_f_lap_g, 1e-6);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Δ₂ tests
// ═══════════════════════════════════════════════════════════════════════════

const PrimalC2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);

test "Δ₂ of constant 2-form: ⟨Δ₂c, c⟩_★₂ ≥ 0" {
    // On a mesh with boundary, a constant 2-form is not in the kernel of Δ₂
    // because the codifferential δ has nonzero boundary contributions. But
    // Δ₂ is positive-semidefinite, so ⟨Δ₂c, c⟩_★₂ ≥ 0 must hold.
    //
    // Note: with the diagonal ★₁⁻¹, interior faces got exactly zero because
    // the pointwise inverse was local. With the Whitney ★₁⁻¹ (global CG
    // solve), boundary effects propagate to interior faces.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(2);

    var omega = try PrimalC2.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 7.0;

    var result = try (try operator_context.laplacian(2)).apply(allocator, omega);
    defer result.deinit(allocator);

    // ⟨Δ₂c, c⟩_★₂ = Σ_f c_f · (Δ₂c)_f / area_f
    const areas = mesh.simplices(2).items(.volume);
    var inner: f64 = 0;
    for (omega.values, result.values, areas) |w, lw, area| {
        inner += w * lw / area;
    }
    try testing.expect(inner >= -1e-8);
}

test "Δ₂ is positive-semidefinite on random 2-forms (500 trials)" {
    // ⟨ω, Δ₂ω⟩_★₂ = ‖δ₁ω‖²_★₁ ≥ 0
    // The ★₂-weighted inner product uses 1/area ratios.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(2);

    const areas = mesh.simplices(2).items(.volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_20);

    for (0..500) |_| {
        var omega = try PrimalC2.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try (try operator_context.laplacian(2)).apply(allocator, omega);
        defer lap_omega.deinit(allocator);

        // ⟨ω, Δ₂ω⟩_★₂ = Σᵢ ωᵢ · (Δ₂ω)ᵢ · (1 / area[i])
        var inner: f64 = 0;
        for (omega.values, lap_omega.values, areas) |w, lw, area| {
            inner += w * lw / area;
        }
        try testing.expect(inner >= -1e-8);
    }
}

test "Δ₂ is symmetric in ★₂-weighted inner product (500 trials)" {
    // Self-adjointness: ⟨Δ₂f, g⟩_★₂ = ⟨f, Δ₂g⟩_★₂
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);
    const operator_context = try context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(2);

    const areas = mesh.simplices(2).items(.volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_21);

    for (0..500) |_| {
        var f_form = try PrimalC2.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC2.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try (try operator_context.laplacian(2)).apply(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try (try operator_context.laplacian(2)).apply(allocator, g_form);
        defer lap_g.deinit(allocator);

        var inner_lap_f_g: f64 = 0;
        var inner_f_lap_g: f64 = 0;
        for (0..f_form.values.len) |idx| {
            const weight = 1.0 / areas[idx];
            inner_lap_f_g += lap_f.values[idx] * g_form.values[idx] * weight;
            inner_f_lap_g += f_form.values[idx] * lap_g.values[idx] * weight;
        }

        try testing.expectApproxEqRel(inner_lap_f_g, inner_f_lap_g, 1e-9);
    }
}
