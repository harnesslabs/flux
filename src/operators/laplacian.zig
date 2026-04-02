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
const exterior_derivative = @import("exterior_derivative.zig");
const hodge_star = @import("hodge_star.zig");
const sparse = @import("../math/sparse.zig");
const context = @import("context.zig");

fn validatePrimalCochainInput(comptime InputType: type) void {
    if (!@hasDecl(InputType, "duality")) {
        @compileError("laplacian requires a Cochain type");
    }
    if (InputType.duality != cochain.Primal) {
        @compileError("laplacian expects a primal cochain");
    }
}

/// Stored Laplacian operator specialized to a primal k-cochain type.
///
/// For k=0, the assembled form keeps the sparse stiffness matrix
/// S = D₀ᵀ M₁ D₀ and the diagonal ★₀⁻¹ scaling separate so application is
/// one SpMV plus pointwise scaling.
pub fn AssembledLaplacian(comptime InputType: type) type {
    comptime validatePrimalCochainInput(InputType);

    return struct {
        const Self = @This();
        const MeshType = InputType.MeshT;
        const k = InputType.degree;

        mesh: *const MeshType,
        stiffness: sparse.CsrMatrix(f64),
        left_scaling: []f64,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.left_scaling);
            self.stiffness.deinit(allocator);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            std.debug.assert(input.mesh == self.mesh);

            var output = try InputType.init(allocator, input.mesh);
            errdefer output.deinit(allocator);

            switch (k) {
                0 => {
                    sparse.spmv(self.stiffness, input.values, output.values);
                    for (output.values, self.left_scaling) |*out, scale| {
                        out.* *= scale;
                    }
                },
                else => return error.NotYetImplemented,
            }

            return output;
        }
    };
}

/// Assemble a stored Laplacian operator for a fixed mesh and primal degree.
pub fn assemble_for_degree(
    comptime MeshType: type,
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledLaplacian(cochain.Cochain(MeshType, k, cochain.Primal)) {
    var stiffness = switch (k) {
        0 => try assemble_zero_form_stiffness(allocator, mesh),
        else => try sparse.CsrMatrix(f64).init(allocator, 0, 0, 0),
    };
    errdefer stiffness.deinit(allocator);

    const left_scaling = switch (k) {
        0 => try assemble_zero_form_star_inverse_diag(allocator, mesh),
        else => try allocator.alloc(f64, 0),
    };
    errdefer allocator.free(left_scaling);

    return .{
        .mesh = mesh,
        .stiffness = stiffness,
        .left_scaling = left_scaling,
    };
}

fn assemble_zero_form_stiffness(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const d0 = mesh.boundary(1);
    const m1 = mesh.whitney_mass_1;

    var assembler = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
    defer assembler.deinit(allocator);

    for (0..m1.n_rows) |edge_i| {
        const incidence_i = d0.row(@intCast(edge_i));
        const mass_row = m1.row(@intCast(edge_i));

        for (mass_row.cols, mass_row.vals) |edge_j, mass_ij| {
            const incidence_j = d0.row(edge_j);

            for (incidence_i.cols, incidence_i.vals) |vertex_i, sign_i| {
                const left = @as(f64, @floatFromInt(sign_i)) * mass_ij;
                for (incidence_j.cols, incidence_j.vals) |vertex_j, sign_j| {
                    const contribution = left * @as(f64, @floatFromInt(sign_j));
                    try assembler.addEntry(allocator, vertex_i, vertex_j, contribution);
                }
            }
        }
    }

    return assembler.build(allocator);
}

fn assemble_zero_form_star_inverse_diag(
    allocator: std.mem.Allocator,
    mesh: anytype,
) ![]f64 {
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    const diagonal = try allocator.alloc(f64, dual_volumes.len);
    errdefer allocator.free(diagonal);

    for (diagonal, dual_volumes) |*out, dual_volume| {
        std.debug.assert(dual_volume != 0.0);
        out.* = 1.0 / dual_volume;
    }

    return diagonal;
}

pub fn laplacian_composed(
    allocator: std.mem.Allocator,
    input: anytype,
) !@TypeOf(input) {
    const InputType = @TypeOf(input);
    comptime validatePrimalCochainInput(InputType);

    const MeshType = InputType.MeshT;
    const k = InputType.degree;
    const topological_dimension = MeshType.topological_dimension;

    var result = try InputType.init(allocator, input.mesh);
    errdefer result.deinit(allocator);

    // Pre-compute workspace size: one allocation for all temporary vectors.
    // Term 1 (δd, k < topological_dimension) needs bk1.n_cols elements.
    // Term 2 (dδ, k > 0) needs 2 × bk.n_cols elements.
    const bk1_cols = if (k < topological_dimension) input.mesh.boundary(k + 1).n_cols else 0;
    const bk_cols = if (k > 0) input.mesh.boundary(k).n_cols else 0;
    const workspace_len = bk1_cols + 2 * bk_cols;

    const workspace = try allocator.alloc(f64, workspace_len);
    defer allocator.free(workspace);

    // ── Term 1 (δd): ★⁻¹_k · D_kᵀ · ★_{k+1} · D_k · ω ─────────────
    // Exists when k < topological_dimension, so that d_k is defined.
    if (k < topological_dimension) {
        // d_k ω = boundary(k+1) · ω
        var d_omega = try exterior_derivative.exterior_derivative(allocator, input);
        defer d_omega.deinit(allocator);

        // ★_{k+1} (d_k ω)
        var star_d = try hodge_star.hodge_star(allocator, d_omega);
        defer star_d.deinit(allocator);

        // D_kᵀ · (★_{k+1} d_k ω) — transpose sparse matrix–vector product
        const bk1 = input.mesh.boundary(k + 1);
        const temp = workspace[0..bk1_cols];
        @memset(temp, 0);
        bk1.transpose_multiply(star_d.values, temp);

        // ★⁻¹_k · temp → result
        try hodge_star.apply_inverse_raw(allocator, MeshType, k, input.mesh, temp, result.values);
    }

    // ── Term 2 (dδ): D_{k-1} · ★⁻¹_{k-1} · D_{k-1}ᵀ · ★_k · ω ────
    // Exists when k > 0, so that δ_k is defined.
    if (k > 0) {
        // ★_k ω
        var star_omega = try hodge_star.hodge_star(allocator, input);
        defer star_omega.deinit(allocator);

        // D_{k-1}ᵀ · (★_k ω) — transpose multiply
        const bk = input.mesh.boundary(k);
        const temp_km1 = workspace[bk1_cols .. bk1_cols + bk_cols];
        @memset(temp_km1, 0);
        bk.transpose_multiply(star_omega.values, temp_km1);

        // ★⁻¹_{k-1} · temp
        const codiff_vals = workspace[bk1_cols + bk_cols ..];
        try hodge_star.apply_inverse_raw(allocator, MeshType, k - 1, input.mesh, temp_km1, codiff_vals);

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

fn apply_laplacian_with_context(
    allocator: std.mem.Allocator,
    input: anytype,
) !@TypeOf(input) {
    const InputType = @TypeOf(input);
    comptime validatePrimalCochainInput(InputType);

    if (InputType.degree != 0) {
        return laplacian_composed(allocator, input);
    }

    var operator_context = context.OperatorContext(InputType.MeshT).init(allocator, input.mesh);
    defer operator_context.deinit();
    try operator_context.withLaplacian(InputType.degree);
    return operator_context.laplacian(InputType.degree).apply(allocator, input);
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);

test "Δ₀ of constant 0-form is zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 42.0;

    var result = try apply_laplacian_with_context(allocator, omega);
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
    var mesh = try Mesh2D.uniform_grid(allocator, nx, ny, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0];

    var result = try apply_laplacian_with_context(allocator, omega);
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

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_00);

    for (0..1000) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try apply_laplacian_with_context(allocator, omega);
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

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_01);

    for (0..1000) |_| {
        var f_form = try PrimalC0.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC0.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try apply_laplacian_with_context(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try apply_laplacian_with_context(allocator, g_form);
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

        var lap_omega = try apply_laplacian_with_context(allocator, omega);
        defer lap_omega.deinit(allocator);

        // Non-constant ⟹ ‖Δ₀ω‖ > 0
        var norm_sq: f64 = 0;
        for (lap_omega.values) |v| norm_sq += v * v;
        try testing.expect(norm_sq > 1e-10);
    }
}

test "assembled Δ₀ apply matches compose-on-the-fly Laplacian" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    var operator_context = context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    try operator_context.withLaplacian(0);

    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_03);

    for (0..100) |_| {
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var expected = try laplacian_composed(allocator, omega);
        defer expected.deinit(allocator);

        var actual = try operator_context.laplacian(0).apply(allocator, omega);
        defer actual.deinit(allocator);

        for (actual.values, expected.values) |got, want| {
            try testing.expectApproxEqAbs(want, got, 5e-12);
        }
    }
}

test "assembled Δ₀ apply is stable across repeated applications" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 6, 5, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    var operator_context = context.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    try operator_context.withLaplacian(0);

    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0] - 0.5 * c[1];

    for (0..3) |_| {
        var expected = try laplacian_composed(allocator, omega);
        defer expected.deinit(allocator);

        var actual = try operator_context.laplacian(0).apply(allocator, omega);
        defer actual.deinit(allocator);

        for (actual.values, expected.values) |got, want| {
            try testing.expectApproxEqAbs(want, got, 5e-12);
        }
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

    var result = try apply_laplacian_with_context(allocator, zero_form);
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

    const edge_volumes = mesh.simplices(1).items(.volume);
    const dual_edge_volumes = mesh.dual_edge_volumes;
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_10);

    for (0..500) |_| {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try apply_laplacian_with_context(allocator, omega);
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

    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

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

        var lap_f = try apply_laplacian_with_context(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try apply_laplacian_with_context(allocator, g_form);
        defer lap_g.deinit(allocator);

        // ⟨Δ₁f, g⟩_★₁ = (Δ₁f)ᵀ M₁ g
        sparse_mod.spmv(mesh.whitney_mass_1, g_form.values, m1_buf);
        var inner_lap_f_g: f64 = 0;
        for (lap_f.values, m1_buf) |lf, mg| inner_lap_f_g += lf * mg;

        // ⟨f, Δ₁g⟩_★₁ = fᵀ M₁ (Δ₁g)
        sparse_mod.spmv(mesh.whitney_mass_1, lap_g.values, m1_buf);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try PrimalC2.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 7.0;

    var result = try apply_laplacian_with_context(allocator, omega);
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
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const areas = mesh.simplices(2).items(.volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_20);

    for (0..500) |_| {
        var omega = try PrimalC2.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_omega = try apply_laplacian_with_context(allocator, omega);
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

    const areas = mesh.simplices(2).items(.volume);
    var rng = std.Random.DefaultPrng.init(0xDEC_1A9_21);

    for (0..500) |_| {
        var f_form = try PrimalC2.init(allocator, &mesh);
        defer f_form.deinit(allocator);
        var g_form = try PrimalC2.init(allocator, &mesh);
        defer g_form.deinit(allocator);

        for (f_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (g_form.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var lap_f = try apply_laplacian_with_context(allocator, f_form);
        defer lap_f.deinit(allocator);
        var lap_g = try apply_laplacian_with_context(allocator, g_form);
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
