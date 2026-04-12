//! Discrete exterior derivative operator.
//!
//! The exterior derivative dₖ: Ωᵏ → Ωᵏ⁺¹ maps k-cochains to (k+1)-cochains.
//! In DEC on a simplicial complex, dₖ is the coboundary operator — the transpose
//! of the boundary operator ∂_{k+1}. On our mesh the boundary matrices are already
//! stored in the coboundary orientation (rows indexed by higher-dimensional cells),
//! so dₖ(ω) is a direct sparse matrix–vector product: boundary_{k+1} · ω.
//!
//! For dual cochains, the dual exterior derivative d̃ₖ uses the transposed primal
//! boundary matrix: d̃ₖ = boundary(n−k)ᵀ. This maps dual k-forms (on primal
//! (n−k)-cells) to dual (k+1)-forms (on primal (n−k−1)-cells).
//!
//! The fundamental identity dd = 0 (boundary_{k+2} · boundary_{k+1} = 0) is an
//! exact algebraic consequence of the oriented incidence structure and holds to
//! machine precision.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");

/// Apply the exterior derivative dₖ to a k-cochain, returning a (k+1)-cochain.
///
/// The degree and mesh type are extracted from the input cochain at comptime.
/// Passing a cochain of the wrong degree to a downstream operator is a compile
/// error because the returned type encodes the degree.
///
/// **Primal cochains** (d on primal Ωᵏ):
///   - d₀: 0-form → 1-form (discrete gradient)
///   - d₁: 1-form → 2-form (discrete curl)
///   - Matrix: boundary(k+1) · ω
///
/// **Dual cochains** (d̃ on dual Ωᵏ):
///   - d̃₀: dual 0-form → dual 1-form
///   - d̃₁: dual 1-form → dual 2-form
///   - Matrix: boundary(n−k)ᵀ · ω
///
/// The dual exterior derivative enables composing ★⁻¹ ∘ d ∘ ★ to express
/// the codifferential δ without manual matrix transposition.
fn ExteriorDerivativeResult(comptime InputType: type) type {
    return cochain.Cochain(InputType.MeshT, InputType.degree + 1, InputType.duality);
}

pub fn AssembledExteriorDerivative(comptime InputType: type) type {
    comptime {
        if (!@hasDecl(InputType, "duality")) {
            @compileError("AssembledExteriorDerivative requires a Cochain type");
        }
    }

    return struct {
        const Self = @This();

        mesh: *const InputType.MeshT,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !ExteriorDerivativeResult(InputType) {
            std.debug.assert(input.mesh == self.mesh);
            return apply_exterior_derivative(allocator, input);
        }
    };
}

pub fn assemble_for_degree(
    comptime MeshType: type,
    comptime Duality: type,
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledExteriorDerivative(cochain.Cochain(MeshType, k, Duality)) {
    _ = allocator;
    return .{ .mesh = mesh };
}

pub fn apply_exterior_derivative(
    allocator: std.mem.Allocator,
    input: anytype,
) !ExteriorDerivativeResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    const k = InputType.degree;
    const topological_dimension = InputType.MeshT.topological_dimension;
    const OutputType = ExteriorDerivativeResult(InputType);

    var output = try OutputType.init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    if (InputType.duality == cochain.Primal) {
        // Primal: dₖ(ω) = boundary(k+1) · ω
        const boundary = input.mesh.boundary(k + 1);
        std.debug.assert(input.values.len == boundary.n_cols);
        boundary.multiply(input.values, output.values);
    } else {
        // Dual: d̃ₖ(ω) = boundary(n−k)ᵀ · ω
        // Dual k-form lives on primal (n−k)-cells; the transpose maps
        // (n−k)-cell values to (n−k−1)-cell values = primal cells for
        // the dual (k+1)-form.
        const boundary = input.mesh.boundary(topological_dimension - k);
        std.debug.assert(input.values.len == boundary.n_rows);

        boundary.transpose_multiply(input.values, output.values);
    }

    return output;
}

pub fn exterior_derivative(
    allocator: std.mem.Allocator,
    input: anytype,
) !ExteriorDerivativeResult(@TypeOf(input)) {
    return apply_exterior_derivative(allocator, input);
}

pub fn assert_closed(
    allocator: std.mem.Allocator,
    input: anytype,
    tolerance: f64,
) !void {
    const InputType = @TypeOf(input);
    comptime {
        if (!@hasDecl(InputType, "degree") or !@hasDecl(InputType, "MeshT")) {
            @compileError("assert_closed requires a Cochain type");
        }
        if (InputType.degree >= InputType.MeshT.topological_dimension) {
            @compileError("assert_closed requires a cochain degree with a defined exterior derivative");
        }
    }

    var derivative = try apply_exterior_derivative(allocator, input);
    defer derivative.deinit(allocator);

    for (derivative.values) |value| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), value, tolerance);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const C0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const C1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const C2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);
const DualC0 = cochain.Cochain(Mesh2D, 0, cochain.Dual);
const DualC1 = cochain.Cochain(Mesh2D, 1, cochain.Dual);
const DualC2 = cochain.Cochain(Mesh2D, 2, cochain.Dual);
const Mesh3D = topology.Mesh(3, 3);
const C3D0 = cochain.Cochain(Mesh3D, 0, cochain.Primal);
const C3D1 = cochain.Cochain(Mesh3D, 1, cochain.Primal);
const C3D2 = cochain.Cochain(Mesh3D, 2, cochain.Primal);
const C3D3 = cochain.Cochain(Mesh3D, 3, cochain.Primal);
const DualC3D0 = cochain.Cochain(Mesh3D, 0, cochain.Dual);
const DualC3D1 = cochain.Cochain(Mesh3D, 1, cochain.Dual);
const DualC3D2 = cochain.Cochain(Mesh3D, 2, cochain.Dual);
const DualC3D3 = cochain.Cochain(Mesh3D, 3, cochain.Dual);

test "d₀ of constant function is zero" {
    // A constant 0-form has zero gradient everywhere.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
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
    var mesh = try Mesh2D.plane(allocator, nx, ny, 3.0, 2.0);
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
    const edge_verts = mesh.simplices(1).items(.vertices);
    for (result.values, edge_verts) |d_val, ev| {
        const expected = coords[ev[1]][0] - coords[ev[0]][0];
        try testing.expectApproxEqAbs(expected, d_val, 1e-14);
    }
}

test "d₁ of d₀ is zero on a specific function" {
    // dd = 0: applying d₁ after d₀ to any 0-form yields the zero 2-form.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
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
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
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
    // There is no d₂ (since topological_dimension = 2), so "dd = 0" for k=1 is vacuously
    // true. Instead, we verify the dual identity: d₁ applied to a random
    // 1-form produces a valid 2-cochain, and that d₁ applied to any exact
    // 1-form (one that is d₀ of something) yields zero — which is the
    // same dd=0 identity tested from the 1-form side.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 6, 5, 3.0, 2.0);
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

test "dd = 0 for random 0-forms on tetrahedral meshes (1000 trials)" {
    const allocator = testing.allocator;
    const mesh_sizes = [_][3]u32{
        .{ 1, 1, 1 },
        .{ 2, 1, 1 },
        .{ 1, 2, 1 },
        .{ 1, 1, 2 },
        .{ 2, 2, 1 },
    };
    var rng = std.Random.DefaultPrng.init(0x3D_DD_00);

    for (mesh_sizes) |size| {
        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, size[0], size[1], size[2], 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        for (0..200) |_| {
            var omega = try C3D0.init(allocator, &mesh);
            defer omega.deinit(allocator);

            for (omega.values) |*value| {
                value.* = rng.random().float(f64) * 200.0 - 100.0;
            }

            var d_omega = try exterior_derivative(allocator, omega);
            defer d_omega.deinit(allocator);

            var dd_omega = try exterior_derivative(allocator, d_omega);
            defer dd_omega.deinit(allocator);

            for (dd_omega.values) |value| {
                try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-11);
            }
        }
    }
}

test "dd = 0 for random 1-forms on tetrahedral meshes (1000 trials)" {
    const allocator = testing.allocator;
    const mesh_sizes = [_][3]u32{
        .{ 1, 1, 1 },
        .{ 2, 1, 1 },
        .{ 1, 2, 1 },
        .{ 1, 1, 2 },
        .{ 2, 2, 1 },
    };
    var rng = std.Random.DefaultPrng.init(0x3D_DD_01);

    for (mesh_sizes) |size| {
        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, size[0], size[1], size[2], 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        for (0..200) |_| {
            var omega = try C3D1.init(allocator, &mesh);
            defer omega.deinit(allocator);

            for (omega.values) |*value| {
                value.* = rng.random().float(f64) * 200.0 - 100.0;
            }

            var d_omega = try exterior_derivative(allocator, omega);
            defer d_omega.deinit(allocator);

            var dd_omega = try exterior_derivative(allocator, d_omega);
            defer dd_omega.deinit(allocator);

            for (dd_omega.values) |value| {
                try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-11);
            }
        }
    }
}

test "d₂B = 0 structural assertion on exact 2-forms in 3D" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0x3D_B_00);

    for (0..1000) |_| {
        var potential = try C3D1.init(allocator, &mesh);
        defer potential.deinit(allocator);

        for (potential.values) |*value| {
            value.* = rng.random().float(f64) * 200.0 - 100.0;
        }

        var magnetic_flux = try exterior_derivative(allocator, potential);
        defer magnetic_flux.deinit(allocator);

        try assert_closed(allocator, magnetic_flux, 1e-11);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Compile-time degree enforcement tests
// ═══════════════════════════════════════════════════════════════════════════

test "compile-time: d₀ returns a 1-cochain, d₁ returns a 2-cochain" {
    comptime {
        // The return type encodes the output degree — passing it where
        // a different degree is expected is a type mismatch.
        const D0_Output = cochain.Cochain(Mesh2D, C0.degree + 1, cochain.Primal);
        const D1_Output = cochain.Cochain(Mesh2D, C1.degree + 1, cochain.Primal);

        try testing.expect(D0_Output == C1);
        try testing.expect(D1_Output == C2);
    }
}

test "compile-time: 3D dₖ returns the next primal degree for k ∈ {0,1,2}" {
    comptime {
        try testing.expect(ExteriorDerivativeResult(C3D0) == C3D1);
        try testing.expect(ExteriorDerivativeResult(C3D1) == C3D2);
        try testing.expect(ExteriorDerivativeResult(C3D2) == C3D3);
    }
}

test "compile-time: Cochain types of different degree are distinct" {
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

// ═══════════════════════════════════════════════════════════════════════════
// Dual exterior derivative tests
// ═══════════════════════════════════════════════════════════════════════════

test "compile-time: d̃₀ maps dual 0-form to dual 1-form" {
    comptime {
        try testing.expect(ExteriorDerivativeResult(DualC0) == DualC1);
    }
}

test "compile-time: d̃₁ maps dual 1-form to dual 2-form" {
    comptime {
        try testing.expect(ExteriorDerivativeResult(DualC1) == DualC2);
    }
}

test "compile-time: 3D d̃ₖ maps dual forms to the next dual degree for k ∈ {0,1,2}" {
    comptime {
        try testing.expect(ExteriorDerivativeResult(DualC3D0) == DualC3D1);
        try testing.expect(ExteriorDerivativeResult(DualC3D1) == DualC3D2);
        try testing.expect(ExteriorDerivativeResult(DualC3D2) == DualC3D3);
    }
}

test "d̃₀ on constant dual 0-form is zero at interior edges" {
    // A constant dual 0-form has d̃₀ = 0 at interior edges (shared by
    // two faces with opposite signs). Boundary edges have only one
    // adjacent face, so d̃₀ is nonzero there — this is expected.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try DualC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 5.0;

    var result = try exterior_derivative(allocator, omega);
    defer result.deinit(allocator);

    // All edges have nonzero dual_length with the barycentric dual.
    // Count nonzero entries at interior edges to verify they are zero.
    const boundary2 = mesh.boundary(2);
    var interior_count: usize = 0;
    for (0..boundary2.n_cols) |edge_idx| {
        // An interior edge appears as a column in boundary(2) with
        // exactly two nonzero entries (two adjacent faces).
        var face_count: usize = 0;
        for (0..boundary2.n_rows) |face_idx| {
            const r = boundary2.row(@intCast(face_idx));
            for (r.cols) |col| {
                if (col == edge_idx) {
                    face_count += 1;
                    break;
                }
            }
        }
        if (face_count == 2) {
            try testing.expectApproxEqAbs(@as(f64, 0), result.values[edge_idx], 1e-15);
            interior_count += 1;
        }
    }
    // Sanity: a 4×3 grid should have many interior edges.
    try testing.expect(interior_count > 0);
}

test "d̃d̃ = 0 for random dual 0-forms (1000 trials)" {
    // The dual exterior derivative satisfies the same dd = 0 identity:
    // d̃₁(d̃₀(ω)) = boundary(1)ᵀ · boundary(2)ᵀ · ω = (boundary(2) · boundary(1))ᵀ · ω = 0
    // because ∂∂ = 0.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_DD_D0);

    for (0..1000) |_| {
        var omega = try DualC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var d_omega = try exterior_derivative(allocator, omega);
        defer d_omega.deinit(allocator);

        var dd_omega = try exterior_derivative(allocator, d_omega);
        defer dd_omega.deinit(allocator);

        for (dd_omega.values) |v| {
            try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-11);
        }
    }
}

test "d̃d̃ = 0 for random dual 0-forms on tetrahedral meshes (1000 trials)" {
    const allocator = testing.allocator;
    const mesh_sizes = [_][3]u32{
        .{ 1, 1, 1 },
        .{ 2, 1, 1 },
        .{ 1, 2, 1 },
        .{ 1, 1, 2 },
        .{ 2, 2, 1 },
    };
    var rng = std.Random.DefaultPrng.init(0x3D_DD_D0);

    for (mesh_sizes) |size| {
        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, size[0], size[1], size[2], 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        for (0..200) |_| {
            var omega = try DualC3D0.init(allocator, &mesh);
            defer omega.deinit(allocator);

            for (omega.values) |*value| {
                value.* = rng.random().float(f64) * 200.0 - 100.0;
            }

            var d_omega = try exterior_derivative(allocator, omega);
            defer d_omega.deinit(allocator);

            var dd_omega = try exterior_derivative(allocator, d_omega);
            defer dd_omega.deinit(allocator);

            for (dd_omega.values) |value| {
                try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-11);
            }
        }
    }
}

test "d̃d̃ = 0 for random dual 1-forms on tetrahedral meshes (1000 trials)" {
    const allocator = testing.allocator;
    const mesh_sizes = [_][3]u32{
        .{ 1, 1, 1 },
        .{ 2, 1, 1 },
        .{ 1, 2, 1 },
        .{ 1, 1, 2 },
        .{ 2, 2, 1 },
    };
    var rng = std.Random.DefaultPrng.init(0x3D_DD_D1);

    for (mesh_sizes) |size| {
        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, size[0], size[1], size[2], 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        for (0..200) |_| {
            var omega = try DualC3D1.init(allocator, &mesh);
            defer omega.deinit(allocator);

            for (omega.values) |*value| {
                value.* = rng.random().float(f64) * 200.0 - 100.0;
            }

            var d_omega = try exterior_derivative(allocator, omega);
            defer d_omega.deinit(allocator);

            var dd_omega = try exterior_derivative(allocator, d_omega);
            defer dd_omega.deinit(allocator);

            for (dd_omega.values) |value| {
                try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-11);
            }
        }
    }
}
