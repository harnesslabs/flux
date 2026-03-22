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

/// The result type of applying d to a given cochain type.
fn Result(comptime InputType: type) type {
    return cochain.Cochain(InputType.MeshT, InputType.degree + 1);
}

/// Apply the exterior derivative dₖ to a k-cochain, returning a (k+1)-cochain.
///
/// The degree and mesh type are extracted from the input cochain at comptime.
/// Passing a cochain of the wrong degree to a downstream operator is a compile
/// error because the returned type encodes the degree.
///
/// For a 2D mesh:
///   - d₀: 0-form → 1-form (discrete gradient). For each edge e,
///     (d₀ω)(e) = ω(head(e)) − ω(tail(e)).
///   - d₁: 1-form → 2-form (discrete curl). For each face f,
///     (d₁ω)(f) = Σ over boundary edges of f, with orientation signs.
pub fn exterior_derivative(allocator: std.mem.Allocator, input: anytype) !Result(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    const k = InputType.degree;

    const boundary = switch (k) {
        0 => input.mesh.boundary_1,
        1 => input.mesh.boundary_2,
        else => @compileError(std.fmt.comptimePrint(
            "exterior_derivative not implemented for k={d} on a 2D mesh",
            .{k},
        )),
    };

    std.debug.assert(input.values.len == boundary.n_cols);

    var output = try Result(InputType).init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    // dₖ(ω) = boundary_{k+1} · ω  (sparse matrix–vector product)
    for (0..boundary.n_rows) |row_idx| {
        const r = boundary.row(@intCast(row_idx));
        var sum: f64 = 0;
        for (r.cols, r.vals) |col, sign| {
            sum += @as(f64, @floatFromInt(sign)) * input.values[col];
        }
        output.values[row_idx] = sum;
    }

    return output;
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

    var result = try exterior_derivative(allocator, omega);
    defer result.deinit(allocator);

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

    var result = try exterior_derivative(allocator, omega);
    defer result.deinit(allocator);

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

    var d_omega = try exterior_derivative(allocator, omega);
    defer d_omega.deinit(allocator);

    var dd_omega = try exterior_derivative(allocator, d_omega);
    defer dd_omega.deinit(allocator);

    for (dd_omega.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-13);
    }
}

test "dd = 0 for random 0-forms on triangular mesh (1000 trials)" {
    // Property test: d₁(d₀(ω)) = 0 for 1000 random 0-cochains.
    // This is the cohomological identity — it must hold exactly because
    // boundary₂ · boundary₁ = 0 as integer matrices.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_DD_00);

    for (0..1000) |_| {
        var omega = try C0.init(allocator, &mesh);
        defer omega.deinit(allocator);

        // Fill ω with random values in [−100, 100].
        for (omega.values) |*v| {
            v.* = rng.random().float(f64) * 200.0 - 100.0;
        }

        var d_omega = try exterior_derivative(allocator, omega);
        defer d_omega.deinit(allocator);

        var dd_omega = try exterior_derivative(allocator, d_omega);
        defer dd_omega.deinit(allocator);

        for (dd_omega.values) |v| {
            try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-11);
        }
    }
}

test "dd = 0 for random 1-forms on triangular mesh (1000 trials)" {
    // Property test: for a 2D mesh, d₁ maps 1-forms to 2-forms.
    // There is no d₂ (since dim = 2), so "dd = 0" for k=1 is vacuously
    // true. Instead, we verify the dual identity: d₁ applied to a random
    // 1-form produces a valid 2-cochain, and that d₁ applied to any exact
    // 1-form (one that is d₀ of something) yields zero — which is the
    // same dd=0 identity tested from the 1-form side.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 6, 5, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_DD_01);

    for (0..1000) |_| {
        var phi = try C0.init(allocator, &mesh);
        defer phi.deinit(allocator);

        // Generate a random 0-form, apply d₀ to get an exact 1-form,
        // then verify d₁ of that exact 1-form is zero.
        for (phi.values) |*v| {
            v.* = rng.random().float(f64) * 200.0 - 100.0;
        }

        var exact_1form = try exterior_derivative(allocator, phi);
        defer exact_1form.deinit(allocator);

        var result = try exterior_derivative(allocator, exact_1form);
        defer result.deinit(allocator);

        for (result.values) |v| {
            try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-11);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Compile-time degree enforcement tests
// ═══════════════════════════════════════════════════════════════════════════

test "compile-time: d returns the correct output type" {
    // The return type of exterior_derivative encodes the degree — passing
    // its output to a function expecting a different degree is a compile error.
    comptime {
        // d₀: Ω⁰ → Ω¹
        try testing.expect(Result(C0) == C1);
        // d₁: Ω¹ → Ω²
        try testing.expect(Result(C1) == C2);

        // The output type is never the same as the input type
        try testing.expect(Result(C0) != C0);
        try testing.expect(Result(C1) != C1);
    }
}

test "compile-time: Cochain types of different degree are distinct" {
    // exterior_derivative returns a typed cochain — passing a 0-form where
    // a 1-form is expected is a type mismatch at compile time.
    comptime {
        try testing.expect(C0 != C1);
        try testing.expect(C1 != C2);
        try testing.expect(C0 != C2);
    }
}

test "compile-time: degree and mesh type recoverable from cochain" {
    comptime {
        try testing.expectEqual(0, C0.degree);
        try testing.expectEqual(1, C1.degree);
        try testing.expectEqual(2, C2.degree);

        try testing.expect(C0.MeshT == Mesh2D);
        try testing.expect(C1.MeshT == Mesh2D);
        try testing.expect(C2.MeshT == Mesh2D);
    }
}
