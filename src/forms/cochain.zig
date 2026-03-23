//! Discrete k-forms (cochains) on simplicial meshes.
//!
//! A k-cochain assigns a real value to every k-cell of a mesh. The degree `k`
//! and mesh type are comptime parameters, so the exterior derivative and other
//! operators can enforce dimensional compatibility as compile errors.

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");

/// Whether a cochain lives on the primal or dual complex.
///
/// The Hodge star ★ maps primal k-forms to dual (n−k)-forms, and ★⁻¹ maps
/// back. Operators that expect a specific duality reject the wrong one at
/// compile time.
pub const Duality = enum { primal, dual };

/// A discrete k-form (cochain) on a simplicial mesh.
///
/// In DEC, a k-cochain assigns a real value to every k-cell of the mesh:
///   - 0-cochain: one value per vertex   (discrete scalar field)
///   - 1-cochain: one value per edge     (discrete circulation/flux)
///   - 2-cochain: one value per face     (discrete flux/density)
///
/// The degree `k`, mesh type, and duality are fixed at comptime, so that
/// operators like the exterior derivative `d` and Hodge star `★` can
/// enforce dimensional and duality compatibility as compile errors.
pub fn Cochain(comptime MeshType: type, comptime k: comptime_int, comptime duality: Duality) type {
    comptime {
        if (!@hasDecl(MeshType, "dimension")) {
            @compileError("Cochain requires a Mesh type with a 'dimension' declaration");
        }
        // A k-cochain on an n-dimensional complex has degree 0 ≤ k ≤ n.
        if (k < 0 or k > MeshType.dimension) {
            @compileError(std.fmt.comptimePrint(
                "Cochain degree k={d} out of range for {d}-dimensional mesh (must be 0 ≤ k ≤ {d})",
                .{ k, MeshType.dimension, MeshType.dimension },
            ));
        }
    }

    return struct {
        const Self = @This();

        /// The mesh type this cochain is defined on.
        pub const MeshT = MeshType;
        /// The degree of this cochain (0 = vertices, 1 = edges, 2 = faces, ...).
        pub const degree = k;
        /// Whether this cochain lives on the primal or dual complex.
        pub const dual = duality;

        /// Coefficient values — one per k-cell.
        values: []f64,
        /// The mesh this cochain is defined on.
        mesh: *const MeshType,

        /// Number of cells this cochain has one value for.
        ///
        /// Primal k-cells map directly: 0 → vertices, 1 → edges, 2 → faces.
        /// Dual k-cells correspond to primal (n−k)-cells: dual 0-cells are
        /// primal faces, dual 1-cells are primal edges, dual 2-cells are
        /// primal vertices.
        pub fn num_cells(mesh: *const MeshType) u32 {
            const effective_degree = if (duality == .dual)
                MeshType.topological_dimension - k
            else
                k;
            return switch (effective_degree) {
                0 => mesh.num_vertices(),
                1 => mesh.num_edges(),
                2 => mesh.num_faces(),
                else => @compileError("unsupported degree for num_cells lookup"),
            };
        }

        /// Allocate a zero-initialized cochain on the given mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            const count = num_cells(mesh);
            const values = try allocator.alloc(f64, count);
            @memset(values, 0);
            return .{ .values = values, .mesh = mesh };
        }

        /// Free the coefficient storage.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.values);
        }

        // ── Arithmetic ──────────────────────────────────────────────────

        /// Pointwise addition: self += other.
        pub fn add(self: *Self, other: Self) void {
            std.debug.assert(self.values.len == other.values.len);
            for (self.values, other.values) |*a, b| {
                a.* += b;
            }
        }

        /// Pointwise subtraction: self -= other.
        pub fn sub(self: *Self, other: Self) void {
            std.debug.assert(self.values.len == other.values.len);
            for (self.values, other.values) |*a, b| {
                a.* -= b;
            }
        }

        /// Scalar multiplication: self *= scalar.
        pub fn scale(self: *Self, scalar: f64) void {
            for (self.values) |*v| {
                v.* *= scalar;
            }
        }

        /// Negate all coefficients in place: self = -self.
        pub fn negate(self: *Self) void {
            for (self.values) |*v| {
                v.* = -v.*;
            }
        }

        /// L² inner product: ⟨self, other⟩ = Σᵢ selfᵢ · otherᵢ.
        /// This is the flat (unweighted) inner product — the Hodge-weighted
        /// version will come with the Hodge star in M2.
        pub fn inner_product(self: Self, other: Self) f64 {
            std.debug.assert(self.values.len == other.values.len);
            var sum: f64 = 0;
            for (self.values, other.values) |a, b| {
                sum += a * b;
            }
            return sum;
        }

        /// Squared L² norm: ‖self‖² = ⟨self, self⟩.
        pub fn norm_squared(self: Self) f64 {
            return self.inner_product(self);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);

test "Cochain(Mesh, 0) has one value per vertex" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_vertices()), c.values.len);
}

test "Cochain(Mesh, 1) has one value per edge" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_edges()), c.values.len);
}

test "Cochain(Mesh, 2) has one value per face" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 2, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_faces()), c.values.len);
}

test "Cochain initializes to zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    for (c.values) |v| {
        try testing.expectEqual(@as(f64, 0), v);
    }
}

test "Cochain degree is accessible at comptime" {
    const C0 = Cochain(Mesh2D, 0, .primal);
    const C1 = Cochain(Mesh2D, 1, .primal);
    const C2 = Cochain(Mesh2D, 2, .primal);

    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 1), C1.degree);
    try testing.expectEqual(@as(comptime_int, 2), C2.degree);
}

test "different degree cochains are distinct types" {
    // This test verifies that Cochain(Mesh, 0) and Cochain(Mesh, 1) are
    // different types — the foundation for compile-time degree enforcement.
    const C0 = Cochain(Mesh2D, 0, .primal);
    const C1 = Cochain(Mesh2D, 1, .primal);
    try testing.expect(C0 != C1);
}

test "primal and dual cochains of the same degree are distinct types" {
    const Primal1 = Cochain(Mesh2D, 1, .primal);
    const Dual1 = Cochain(Mesh2D, 1, .dual);
    try testing.expect(Primal1 != Dual1);

    // Duality is recoverable at comptime.
    try testing.expect(Primal1.dual == .primal);
    try testing.expect(Dual1.dual == .dual);
}

test "Cochain degree bounded by mesh dimension" {
    // Mesh(3) supports 0-, 1-, 2-, and 3-cochains.
    const Mesh3D = topology.Mesh(3);
    const C0 = Cochain(Mesh3D, 0, .primal);
    const C3 = Cochain(Mesh3D, 3, .primal);
    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 3), C3.degree);
}

// ═══════════════════════════════════════════════════════════════════════════
// Arithmetic tests
// ═══════════════════════════════════════════════════════════════════════════

test "add two 0-cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var a = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer b.deinit(allocator);

    for (a.values, 0..) |*v, i| v.* = @floatFromInt(i);
    for (b.values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) * 10.0;

    a.add(b);

    for (a.values, 0..) |v, i| {
        const fi: f64 = @floatFromInt(i);
        try testing.expectApproxEqAbs(fi + fi * 10.0, v, 1e-15);
    }
}

test "sub two 1-cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var a = try Cochain(Mesh2D, 1, .primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 1, .primal).init(allocator, &mesh);
    defer b.deinit(allocator);

    for (a.values, 0..) |*v, i| v.* = @floatFromInt(i);
    for (b.values) |*v| v.* = 3.0;

    a.sub(b);

    for (a.values, 0..) |v, i| {
        const fi: f64 = @floatFromInt(i);
        try testing.expectApproxEqAbs(fi - 3.0, v, 1e-15);
    }
}

test "scale a cochain" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    for (c.values, 0..) |*v, i| v.* = @floatFromInt(i);
    c.scale(2.5);

    for (c.values, 0..) |v, i| {
        const fi: f64 = @floatFromInt(i);
        try testing.expectApproxEqAbs(fi * 2.5, v, 1e-15);
    }
}

test "negate a cochain" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, .primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    for (c.values, 0..) |*v, i| v.* = @floatFromInt(i);
    c.negate();

    for (c.values, 0..) |v, i| {
        const fi: f64 = @floatFromInt(i);
        try testing.expectApproxEqAbs(-fi, v, 1e-15);
    }
}

test "inner product and norm" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // 1×1 grid has 4 vertices
    var a = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 0, .primal).init(allocator, &mesh);
    defer b.deinit(allocator);

    a.values[0] = 1.0;
    a.values[1] = 2.0;
    a.values[2] = 3.0;
    a.values[3] = 4.0;
    b.values[0] = 2.0;
    b.values[1] = 0.0;
    b.values[2] = -1.0;
    b.values[3] = 3.0;

    // ⟨a, b⟩ = 1·2 + 2·0 + 3·(−1) + 4·3 = 2 + 0 − 3 + 12 = 11
    try testing.expectApproxEqAbs(@as(f64, 11.0), a.inner_product(b), 1e-15);

    // ‖a‖² = 1 + 4 + 9 + 16 = 30
    try testing.expectApproxEqAbs(@as(f64, 30.0), a.norm_squared(), 1e-15);
}
