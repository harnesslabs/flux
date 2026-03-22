//! Discrete exterior derivative operator.
//!
//! The exterior derivative dₖ: Ωᵏ → Ωᵏ⁺¹ maps k-cochains to (k+1)-cochains.
//! In DEC on a simplicial complex, dₖ is the coboundary operator — the transpose
//! of the boundary operator ∂_{k+1}. On our mesh the boundary matrices are already
//! stored in the coboundary orientation (rows indexed by higher-dimensional cells),
//! so dₖ(ω) is a direct sparse matrix–vector product: boundary_{k+1} · ω.
//!
//! The fundamental identity dd = 0 (boundary_{k+2} · boundary_{k+1} = 0) is an
//! exact algebraic consequence of the oriented incidence structure and holds to
//! machine precision.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");

/// Apply the exterior derivative dₖ to a k-cochain, producing a (k+1)-cochain.
///
/// For a 2D mesh:
///   - d₀: 0-form → 1-form (discrete gradient). For each edge e,
///     (d₀ω)(e) = ω(head(e)) − ω(tail(e)).
///   - d₁: 1-form → 2-form (discrete curl). For each face f,
///     (d₁ω)(f) = Σ over boundary edges of f, with orientation signs.
///
/// The output cochain must be pre-allocated with the correct number of
/// (k+1)-cells. The caller owns both input and output memory.
pub fn exterior_derivative(
    comptime MeshType: type,
    comptime k: comptime_int,
    input: cochain.Cochain(MeshType, k),
    output: *cochain.Cochain(MeshType, k + 1),
) void {
    // Degree bounds are enforced at comptime by Cochain — if k+1 > dimension,
    // the Cochain(MeshType, k+1) parameter itself will not compile.

    const boundary = switch (k) {
        0 => input.mesh.boundary_1,
        1 => input.mesh.boundary_2,
        else => @compileError(std.fmt.comptimePrint(
            "exterior_derivative not implemented for k={d} on a 2D mesh",
            .{k},
        )),
    };

    std.debug.assert(output.values.len == boundary.n_rows);
    std.debug.assert(input.values.len == boundary.n_cols);

    // dₖ(ω) = boundary_{k+1} · ω  (sparse matrix–vector product)
    for (0..boundary.n_rows) |row_idx| {
        const r = boundary.row(@intCast(row_idx));
        var sum: f64 = 0;
        for (r.cols, r.vals) |col, sign| {
            sum += @as(f64, @floatFromInt(sign)) * input.values[col];
        }
        output.values[row_idx] = sum;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);
const C0 = cochain.Cochain(Mesh2D, 0);
const C1 = cochain.Cochain(Mesh2D, 1);
const C2 = cochain.Cochain(Mesh2D, 2);

test "d₀ of constant function is zero" {
    // A constant 0-form has zero gradient everywhere.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try C0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 7.0; // constant function

    var result = try C1.init(allocator, &mesh);
    defer result.deinit(allocator);

    exterior_derivative(Mesh2D, 0, omega, &result);

    for (result.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-15);
    }
}

test "d₀ on linear function f(x,y) = x" {
    // On a uniform grid, d₀ of f(x,y) = x should give the x-component
    // of the gradient integrated along each edge.
    const allocator = testing.allocator;
    const nx: u32 = 3;
    const ny: u32 = 2;
    var mesh = try Mesh2D.uniform_grid(allocator, nx, ny, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var omega = try C0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    // Set ω(v) = x-coordinate of vertex v
    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*val, c| {
        val.* = c[0];
    }

    var result = try C1.init(allocator, &mesh);
    defer result.deinit(allocator);

    exterior_derivative(Mesh2D, 0, omega, &result);

    // Verify: (d₀ω)(e) = x(head) − x(tail) for each edge
    const edge_verts = mesh.edges.slice().items(.vertices);
    for (result.values, edge_verts) |d_val, ev| {
        const expected = coords[ev[1]][0] - coords[ev[0]][0];
        try testing.expectApproxEqAbs(expected, d_val, 1e-14);
    }
}

test "d₁ of d₀ is zero on a specific function" {
    // dd = 0: applying d₁ after d₀ to any 0-form yields the zero 2-form.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var omega = try C0.init(allocator, &mesh);
    defer omega.deinit(allocator);

    // Use f(x,y) = x² + 3xy − y (nonlinear, nontrivial)
    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*val, c| {
        val.* = c[0] * c[0] + 3.0 * c[0] * c[1] - c[1];
    }

    var d_omega = try C1.init(allocator, &mesh);
    defer d_omega.deinit(allocator);
    exterior_derivative(Mesh2D, 0, omega, &d_omega);

    var dd_omega = try C2.init(allocator, &mesh);
    defer dd_omega.deinit(allocator);
    exterior_derivative(Mesh2D, 1, d_omega, &dd_omega);

    for (dd_omega.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-13);
    }
}
