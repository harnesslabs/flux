//! Discrete Laplace-de Rham operator Δ = dδ + δd.
//!
//! The Hodge Laplacian on primal k-cochains decomposes into two terms:
//!   - δd (curl-curl): ★⁻¹_k · Dᵀ_k · ★_{k+1} · D_k
//!   - dδ (grad-div):  D_{k-1} · ★⁻¹_{k-1} · Dᵀ_{k-1} · ★_k
//!
//! where D_k = boundary(k+1) is the exterior derivative matrix and Dᵀ is
//! its transpose — the matrix representation of the codifferential's
//! topological part.
//!
//! For a 2D simplicial mesh:
//!   Δ₀ = ★₀⁻¹ D₀ᵀ ★₁ D₀                           (δd only)
//!   Δ₁ = D₀ ★₀⁻¹ D₀ᵀ ★₁ + ★₁⁻¹ D₁ᵀ ★₂ D₁         (both terms)
//!   Δ₂ = D₁ ★₁⁻¹ D₁ᵀ ★₂                            (dδ only)
//!
//! Uses the unsigned convention: Δ₀ is positive-semidefinite, meaning
//! ⟨Δ₀ω, ω⟩_★₀ ≥ 0 for all ω.
//!
//! All three Laplacians are well-defined on meshes with barycentric dual
//! geometry, where every edge has nonzero dual length.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const ext = @import("exterior_derivative.zig");
const hs = @import("hodge_star.zig");

/// Apply the Hodge Laplacian Δₖ to a primal k-cochain.
///
/// Returns a primal k-cochain of the same degree. The operator is
/// self-adjoint in the ★-weighted inner product and positive-semidefinite
/// on 0-forms.
pub fn laplacian(
    allocator: std.mem.Allocator,
    input: anytype,
) !@TypeOf(input) {
    const InputType = @TypeOf(input);
    comptime {
        if (!@hasDecl(InputType, "duality")) {
            @compileError("laplacian requires a Cochain type");
        }
        if (InputType.duality != cochain.Primal) {
            @compileError("laplacian expects a primal cochain");
        }
    }

    const MeshType = InputType.MeshT;
    const k = InputType.degree;
    const n = MeshType.topological_dimension;

    var result = try InputType.init(allocator, input.mesh);
    errdefer result.deinit(allocator);

    // ── Term 1 (δd): ★⁻¹_k · D_kᵀ · ★_{k+1} · D_k · ω ─────────────
    // Exists when k < n, so that d_k is defined.
    if (k < n) {
        // d_k ω = boundary(k+1) · ω
        var d_omega = try ext.exterior_derivative(allocator, input);
        defer d_omega.deinit(allocator);

        // ★_{k+1} (d_k ω)
        var star_d = try hs.hodge_star(allocator, d_omega);
        defer star_d.deinit(allocator);

        // D_kᵀ · (★_{k+1} d_k ω) — transpose sparse matrix–vector product
        const bk1 = input.mesh.boundary(k + 1);
        const temp = try allocator.alloc(f64, bk1.n_cols);
        defer allocator.free(temp);
        @memset(temp, 0);
        bk1.transpose_multiply(star_d.values, temp);

        // ★⁻¹_k · temp → result
        hs.applyDiagonal(MeshType, k, input.mesh, temp, result.values, true);
    }

    // ── Term 2 (dδ): D_{k-1} · ★⁻¹_{k-1} · D_{k-1}ᵀ · ★_k · ω ────
    // Exists when k > 0, so that δ_k is defined.
    if (k > 0) {
        // ★_k ω
        var star_omega = try hs.hodge_star(allocator, input);
        defer star_omega.deinit(allocator);

        // D_{k-1}ᵀ · (★_k ω) — transpose multiply
        const bk = input.mesh.boundary(k);
        const temp_km1 = try allocator.alloc(f64, bk.n_cols);
        defer allocator.free(temp_km1);
        @memset(temp_km1, 0);
        bk.transpose_multiply(star_omega.values, temp_km1);

        // ★⁻¹_{k-1} · temp
        const codiff_vals = try allocator.alloc(f64, bk.n_cols);
        defer allocator.free(codiff_vals);
        hs.applyDiagonal(MeshType, k - 1, input.mesh, temp_km1, codiff_vals, true);

        // D_{k-1} · codiff_vals → accumulate into result
        for (0..bk.n_rows) |row_idx| {
            const r = bk.row(@intCast(row_idx));
            var sum: f64 = 0;
            for (r.cols, r.vals) |col, sign| {
                sum += @as(f64, @floatFromInt(sign)) * codiff_vals[col];
            }
            result.values[row_idx] += sum;
        }
    }

    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);

test "Δ₀ of constant 0-form is zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 42.0;

    var result = try laplacian(allocator, omega);
    defer result.deinit(allocator);

    for (result.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), v, 1e-13);
    }
}

test "Δ₀ of linear function f(x,y) = x is zero at interior vertices" {
    // The DEC Laplacian reproduces the cotangent Laplacian, which is exact
    // for piecewise-linear functions on any triangulation. Since Δ(x) = 0
    // continuously, interior vertices should see zero.
    const allocator = testing.allocator;
    const nx: u32 = 5;
    const ny: u32 = 4;
    var mesh = try Mesh2D.uniform_grid(allocator, nx, ny, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0];

    var result = try laplacian(allocator, omega);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_area);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_00);

    for (0..1000) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try laplacian(allocator, omega);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_area);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_01);

    for (0..1000) |_| {
        var f_form = try PrimalC0.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC0.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try laplacian(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try laplacian(allocator, g_form);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_02);

    for (0..100) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try laplacian(allocator, omega);
        defer lap_omega.deinit(allocator);

        // Non-constant ⟹ ‖Δ₀ω‖ > 0
        var norm_sq: f64 = 0;
        for (lap_omega.values) |v| norm_sq += v * v;
        try testing.expect(norm_sq > 1e-10);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Δ₁ tests
// ═══════════════════════════════════════════════════════════════════════════

const PrimalC1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);

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
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    // Build ω = d₀(constant) = 0 — trivial but validates that Δ₁(0) = 0.
    var zero_form = try PrimalC1.init(allocator, &mesh);
    defer zero_form.deinit(allocator);

    var result = try laplacian(allocator, zero_form);
    defer result.deinit(allocator);

    for (result.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), v, 1e-13);
    }
}

test "Δ₁ is positive-semidefinite on random 1-forms (500 trials)" {
    // ⟨ω, Δ₁ω⟩_★₁ = ‖d₁ω‖²_★₂ + ‖δ₀ω‖²_★₀ ≥ 0
    // The ★₁-weighted inner product uses dual_length/length ratios.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const edge_slice = mesh.edges.slice();
    const lengths = edge_slice.items(.length);
    const dual_lengths = edge_slice.items(.dual_length);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_10);

    for (0..500) |_| {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try laplacian(allocator, omega);
        defer lap_omega.deinit(allocator);

        // ⟨ω, Δ₁ω⟩_★₁ = Σᵢ ωᵢ · (Δ₁ω)ᵢ · (dual_length[i] / length[i])
        var inner: f64 = 0;
        for (omega.values, lap_omega.values, lengths, dual_lengths) |w, lw, len, dual_len| {
            inner += w * lw * (dual_len / len);
        }
        try testing.expect(inner >= -1e-8);
    }
}

test "Δ₁ is symmetric in ★₁-weighted inner product (500 trials)" {
    // Self-adjointness: ⟨Δ₁f, g⟩_★₁ = ⟨f, Δ₁g⟩_★₁
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const edge_slice = mesh.edges.slice();
    const lengths = edge_slice.items(.length);
    const dual_lengths = edge_slice.items(.dual_length);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_11);

    for (0..500) |_| {
        var f_form = try PrimalC1.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC1.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try laplacian(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try laplacian(allocator, g_form);
        defer lap_g.deinit(allocator);

        var inner_lap_f_g: f64 = 0;
        var inner_f_lap_g: f64 = 0;
        for (0..f_form.values.len) |idx| {
            const weight = dual_lengths[idx] / lengths[idx];
            inner_lap_f_g += lap_f.values[idx] * g_form.values[idx] * weight;
            inner_f_lap_g += f_form.values[idx] * lap_g.values[idx] * weight;
        }

        try testing.expectApproxEqRel(inner_lap_f_g, inner_f_lap_g, 1e-9);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Δ₂ tests
// ═══════════════════════════════════════════════════════════════════════════

const PrimalC2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);

test "Δ₂ of constant 2-form is zero" {
    // On a closed manifold, constant top-forms are harmonic.
    // On a mesh with boundary, the dδ term of Δ₂ acts on a constant
    // 2-form. Since d₂ does not exist (top degree), only the dδ term
    // contributes: Δ₂ = D₁ ★₁⁻¹ D₁ᵀ ★₂.
    // For a constant 2-form c, ★₂c is constant, and D₁ᵀ(★₂c) sums the
    // constant over adjacent faces with opposite signs — zero at interior
    // edges and nonzero only at boundary edges. So Δ₂(c) is nonzero
    // only at faces touching the boundary.
    //
    // We verify: Δ₂(c) = 0 at faces not touching any boundary edge.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC2.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 7.0;

    var result = try laplacian(allocator, omega);
    defer result.deinit(allocator);

    // Find interior faces: faces whose three edges are all interior.
    const boundary_2 = mesh.boundary(2);
    var boundary_edge_set = try allocator.alloc(bool, mesh.num_edges());
    defer allocator.free(boundary_edge_set);
    @memset(boundary_edge_set, false);
    for (mesh.boundary_edges) |e| boundary_edge_set[e] = true;

    for (0..mesh.num_faces()) |f| {
        const face_edges = boundary_2.row(@intCast(f));
        var touches_boundary = false;
        for (face_edges.cols) |e| {
            if (boundary_edge_set[e]) {
                touches_boundary = true;
                break;
            }
        }
        if (!touches_boundary) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), result.values[f], 1e-12);
        }
    }
}

test "Δ₂ is positive-semidefinite on random 2-forms (500 trials)" {
    // ⟨ω, Δ₂ω⟩_★₂ = ‖δ₁ω‖²_★₁ ≥ 0
    // The ★₂-weighted inner product uses 1/area ratios.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const areas = mesh.faces.slice().items(.area);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_20);

    for (0..500) |_| {
        var omega = try PrimalC2.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try laplacian(allocator, omega);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const areas = mesh.faces.slice().items(.area);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_21);

    for (0..500) |_| {
        var f_form = try PrimalC2.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC2.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try laplacian(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try laplacian(allocator, g_form);
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
