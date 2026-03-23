//! Discrete Hodge star operator and its inverse.
//!
//! The diagonal Hodge star ★ₖ maps primal k-cochains to dual (n−k)-cochains
//! via the circumcentric dual mesh. For a 2D simplicial mesh:
//!
//!   (★₀)ᵢᵢ = |dual cell of vertex i|  = dual_area[i]
//!   (★₁)ᵢᵢ = |dual edge i| / |edge i| = dual_length[i] / length[i]
//!   (★₂)ᵢᵢ = 1 / |face i|             = 1 / area[i]
//!
//! The inverse ★⁻¹ is the element-wise reciprocal of the diagonal, mapping
//! dual (n−k)-cochains back to primal k-cochains.
//!
//! **Known degeneracy:** On uniform grids with SW→NE diagonal triangulation,
//! both triangles sharing a diagonal have the same circumcenter, giving
//! dual_length = 0 for diagonal edges. The Hodge star maps those entries to
//! zero; the inverse is undefined there and will panic.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");

/// Apply the Hodge star ★ₖ to a primal k-cochain, returning a dual (n−k)-cochain.
///
/// The Hodge star is a diagonal operator: each output value is the input value
/// scaled by the ratio |dual cell| / |primal cell| for the corresponding mesh
/// entity. This ratio encodes the local metric geometry of the mesh.
pub fn hodge_star(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInput(InputType);

    const k = InputType.degree;
    const n = InputType.MeshT.topological_dimension;
    const OutputType = cochain.Cochain(InputType.MeshT, n - k, cochain.Dual);

    var output = try OutputType.init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    applyDiagonal(InputType.MeshT, k, input.mesh, input.values, output.values, false);

    return output;
}

/// Apply the inverse Hodge star ★⁻¹ to a dual (n−k)-cochain, returning a
/// primal k-cochain.
///
/// The inverse is the element-wise reciprocal of the Hodge star diagonal.
/// Panics (via assertion) if any diagonal entry is zero — this indicates a
/// degenerate mesh element (e.g., a diagonal edge on a uniform grid whose
/// dual length is zero).
pub fn hodge_star_inverse(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInverseInput(InputType);

    const dual_degree = InputType.degree;
    const n = InputType.MeshT.topological_dimension;
    const primal_degree = n - dual_degree;
    const OutputType = cochain.Cochain(InputType.MeshT, primal_degree, cochain.Primal);

    var output = try OutputType.init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    applyDiagonal(InputType.MeshT, primal_degree, input.mesh, input.values, output.values, true);

    return output;
}

// ── Return type helpers ──────────────────────────────────────────────────

fn HodgeStarResult(comptime InputType: type) type {
    const n = InputType.MeshT.topological_dimension;
    return cochain.Cochain(InputType.MeshT, n - InputType.degree, cochain.Dual);
}

fn HodgeStarInverseResult(comptime InputType: type) type {
    const n = InputType.MeshT.topological_dimension;
    return cochain.Cochain(InputType.MeshT, n - InputType.degree, cochain.Primal);
}

// ── Comptime validation ──────────────────────────────────────────────────

fn validateHodgeStarInput(comptime T: type) void {
    if (!@hasDecl(T, "duality")) {
        @compileError("hodge_star requires a Cochain type");
    }
    if (T.duality != cochain.Primal) {
        @compileError("hodge_star expects a primal cochain (use hodge_star_inverse for dual → primal)");
    }
}

fn validateHodgeStarInverseInput(comptime T: type) void {
    if (!@hasDecl(T, "duality")) {
        @compileError("hodge_star_inverse requires a Cochain type");
    }
    if (T.duality != cochain.Dual) {
        @compileError("hodge_star_inverse expects a dual cochain (use hodge_star for primal → dual)");
    }
}

// ── Diagonal application ─────────────────────────────────────────────────

/// Apply the Hodge star diagonal (or its inverse) for a given primal degree.
///
/// When `invert` is false: output[i] = ratio[i] * input[i]
/// When `invert` is true:  output[i] = input[i] / ratio[i]  (asserts ratio ≠ 0)
fn applyDiagonal(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
    comptime invert: bool,
) void {
    switch (primal_degree) {
        // ★₀: ratio = dual_area[i]  (0-cell has unit "volume")
        0 => {
            const dual_areas = mesh.vertices.slice().items(.dual_area);
            for (output, input, dual_areas) |*out, in_val, ratio| {
                if (invert) {
                    std.debug.assert(ratio != 0.0);
                    out.* = in_val / ratio;
                } else {
                    out.* = ratio * in_val;
                }
            }
        },
        // ★₁: ratio = dual_length[i] / length[i]
        1 => {
            const edge_slice = mesh.edges.slice();
            const lengths = edge_slice.items(.length);
            const dual_lengths = edge_slice.items(.dual_length);
            for (output, input, lengths, dual_lengths) |*out, in_val, len, dual_len| {
                if (invert) {
                    std.debug.assert(dual_len != 0.0);
                    out.* = (len / dual_len) * in_val;
                } else {
                    out.* = (dual_len / len) * in_val;
                }
            }
        },
        // ★₂: ratio = 1 / area[i]
        2 => {
            const areas = mesh.faces.slice().items(.area);
            for (output, input, areas) |*out, in_val, area| {
                if (invert) {
                    out.* = area * in_val;
                } else {
                    std.debug.assert(area != 0.0);
                    out.* = in_val / area;
                }
            }
        },
        else => @compileError("unsupported degree for Hodge star diagonal"),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const PrimalC1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const PrimalC2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);
const DualC0 = cochain.Cochain(Mesh2D, 0, cochain.Dual);
const DualC1 = cochain.Cochain(Mesh2D, 1, cochain.Dual);
const DualC2 = cochain.Cochain(Mesh2D, 2, cochain.Dual);

// ── Compile-time type checks ─────────────────────────────────────────────

test "compile-time: ★₀ maps primal 0-form to dual 2-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC0) == DualC2);
    }
}

test "compile-time: ★₁ maps primal 1-form to dual 1-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC1) == DualC1);
    }
}

test "compile-time: ★₂ maps primal 2-form to dual 0-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC2) == DualC0);
    }
}

test "compile-time: ★⁻¹ maps dual (n−k)-form back to primal k-form" {
    comptime {
        try testing.expect(HodgeStarInverseResult(DualC2) == PrimalC0);
        try testing.expect(HodgeStarInverseResult(DualC1) == PrimalC1);
        try testing.expect(HodgeStarInverseResult(DualC0) == PrimalC2);
    }
}

// ── Correctness: Hodge star values ───────────────────────────────────────

test "★₀ scales by dual vertex area" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_area);
    for (result.values, dual_areas) |r, expected| {
        try testing.expectApproxEqAbs(expected, r, 1e-15);
    }
}

test "★₁ scales by dual_length / length" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC1.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    const edge_slice = mesh.edges.slice();
    const lengths = edge_slice.items(.length);
    const dual_lengths = edge_slice.items(.dual_length);
    for (result.values, lengths, dual_lengths) |r, len, dual_len| {
        try testing.expectApproxEqAbs(dual_len / len, r, 1e-15);
    }
}

test "★₂ scales by 1 / face area" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 2, 2.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC2.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    const areas = mesh.faces.slice().items(.area);
    for (result.values, areas) |r, area| {
        try testing.expectApproxEqAbs(1.0 / area, r, 1e-15);
    }
}

// ── Invariant: ★⁻¹ ∘ ★ = identity ───────────────────────────────────────

test "Hodge star inverse is exact for all degrees" {
    // For 1000 random primal cochains of each degree, ★⁻¹(★(ω)) = ω
    // must hold to machine precision.
    //
    // k=0 and k=2 are tested via the full round-trip (all entries non-degenerate).
    // k=1 is tested entry-by-entry, skipping degenerate diagonal edges where
    // dual_length = 0 on the uniform grid.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    // The round-trip involves two floating-point operations (multiply + divide),
    // so we use relative tolerance. Values are in [−100, 100]; 1e-14 is well
    // within the ~15 significant digits of f64.
    const tolerance = 1e-14;

    // ── k = 0: full round-trip ───────────────────────────────────────
    {
        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_00);
        for (0..1000) |_| {
            var omega = try PrimalC0.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, tolerance);
            }
        }
    }

    // ── k = 1: entry-by-entry on non-degenerate edges ────────────────
    {
        const edge_slice = mesh.edges.slice();
        const lengths = edge_slice.items(.length);
        const dual_lengths = edge_slice.items(.dual_length);

        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_01);
        for (0..1000) |_| {
            var omega = try PrimalC1.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            for (omega.values, starred.values, lengths, dual_lengths) |original, star_val, len, dual_len| {
                if (dual_len == 0.0) {
                    // Degenerate: ★₁ maps to zero, inverse is undefined.
                    try testing.expectApproxEqAbs(@as(f64, 0.0), star_val, tolerance);
                } else {
                    // Non-degenerate: manual inverse recovers the original.
                    const recovered = star_val * (len / dual_len);
                    try testing.expectApproxEqRel(original, recovered, tolerance);
                }
            }
        }
    }

    // ── k = 2: full round-trip ───────────────────────────────────────
    {
        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_02);
        for (0..1000) |_| {
            var omega = try PrimalC2.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, tolerance);
            }
        }
    }
}
