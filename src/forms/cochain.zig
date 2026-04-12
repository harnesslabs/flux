//! Discrete k-forms (cochains) on simplicial meshes.
//!
//! A k-cochain assigns a real value to every k-cell of a mesh. The degree `k`
//! and mesh type are comptime parameters, so the exterior derivative and other
//! operators can enforce dimensional compatibility as compile errors.

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");

/// Marker type: cochain lives on the primal complex.
pub const Primal = struct {};

/// Marker type: cochain lives on the dual complex.
/// The Hodge star ★ maps primal k-forms to dual (n−k)-forms, and ★⁻¹ maps
/// back. Operators that expect a specific duality reject the wrong one at
/// compile time via these marker types.
pub const Dual = struct {};

fn addAssignScalar(lhs: []f64, rhs: []const f64) void {
    std.debug.assert(lhs.len == rhs.len);
    for (lhs, rhs) |*left, right| {
        left.* += right;
    }
}

fn addAssignSimd(lhs: []f64, rhs: []const f64) void {
    std.debug.assert(lhs.len == rhs.len);
    if (std.simd.suggestVectorLength(f64)) |lane_count| {
        const Block = @Vector(lane_count, f64);
        const simd_len = lhs.len - (lhs.len % lane_count);

        var i: usize = 0;
        while (i < simd_len) : (i += lane_count) {
            const left_block: Block = lhs[i..][0..lane_count].*;
            const right_block: Block = rhs[i..][0..lane_count].*;
            const sum_array: [lane_count]f64 = left_block + right_block;
            lhs[i..][0..lane_count].* = sum_array;
        }

        addAssignScalar(lhs[simd_len..], rhs[simd_len..]);
        return;
    }

    addAssignScalar(lhs, rhs);
}

fn scaleInPlaceScalar(values: []f64, scalar: f64) void {
    for (values) |*value| {
        value.* *= scalar;
    }
}

fn scaleInPlaceSimd(values: []f64, scalar: f64) void {
    if (std.simd.suggestVectorLength(f64)) |lane_count| {
        const Block = @Vector(lane_count, f64);
        const factor: Block = @splat(scalar);
        const simd_len = values.len - (values.len % lane_count);

        var i: usize = 0;
        while (i < simd_len) : (i += lane_count) {
            const block: Block = values[i..][0..lane_count].*;
            const scaled_array: [lane_count]f64 = block * factor;
            values[i..][0..lane_count].* = scaled_array;
        }

        scaleInPlaceScalar(values[simd_len..], scalar);
        return;
    }

    scaleInPlaceScalar(values, scalar);
}

fn negateInPlaceScalar(values: []f64) void {
    for (values) |*value| {
        value.* = -value.*;
    }
}

fn negateInPlaceSimd(values: []f64) void {
    if (std.simd.suggestVectorLength(f64)) |lane_count| {
        const Block = @Vector(lane_count, f64);
        const simd_len = values.len - (values.len % lane_count);

        var i: usize = 0;
        while (i < simd_len) : (i += lane_count) {
            const block: Block = values[i..][0..lane_count].*;
            const negated_array: [lane_count]f64 = -block;
            values[i..][0..lane_count].* = negated_array;
        }

        negateInPlaceScalar(values[simd_len..]);
        return;
    }

    negateInPlaceScalar(values);
}

fn innerProductScalar(lhs: []const f64, rhs: []const f64) f64 {
    std.debug.assert(lhs.len == rhs.len);
    var sum: f64 = 0;
    for (lhs, rhs) |left, right| {
        sum += left * right;
    }
    return sum;
}

fn innerProductSimd(lhs: []const f64, rhs: []const f64) f64 {
    std.debug.assert(lhs.len == rhs.len);
    if (std.simd.suggestVectorLength(f64)) |lane_count| {
        const Block = @Vector(lane_count, f64);
        const simd_len = lhs.len - (lhs.len % lane_count);

        var sum_block: Block = @splat(0.0);
        var i: usize = 0;
        while (i < simd_len) : (i += lane_count) {
            const left_block: Block = lhs[i..][0..lane_count].*;
            const right_block: Block = rhs[i..][0..lane_count].*;
            sum_block += left_block * right_block;
        }

        return @reduce(.Add, sum_block) + innerProductScalar(lhs[simd_len..], rhs[simd_len..]);
    }

    return innerProductScalar(lhs, rhs);
}

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
pub fn Cochain(comptime MeshType: type, comptime k: comptime_int, comptime D: type) type {
    comptime {
        if (D != Primal and D != Dual) {
            @compileError("Cochain duality parameter must be Primal or Dual");
        }
        if (!@hasDecl(MeshType, "topological_dimension")) {
            @compileError("Cochain requires a Mesh type with a 'topological_dimension' declaration");
        }
        // A k-cochain on a topological_dimension-dimensional complex has
        // degree 0 ≤ k ≤ topological_dimension.
        if (k < 0 or k > MeshType.topological_dimension) {
            @compileError(std.fmt.comptimePrint(
                "Cochain degree k={d} out of range for {d}-dimensional mesh (must be 0 ≤ k ≤ {d})",
                .{ k, MeshType.topological_dimension, MeshType.topological_dimension },
            ));
        }
    }

    return struct {
        const Self = @This();

        /// The mesh type this cochain is defined on.
        pub const MeshT = MeshType;
        /// The degree of this cochain (0 = vertices, 1 = edges, 2 = faces, ...).
        pub const degree = k;
        /// Whether this cochain lives on the primal or dual complex (Primal or Dual).
        pub const duality = D;

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
            const effective_degree = if (D == Dual)
                MeshType.topological_dimension - k
            else
                k;
            return mesh.num_cells(effective_degree);
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
            addAssignSimd(self.values, other.values);
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
            scaleInPlaceSimd(self.values, scalar);
        }

        /// Negate all coefficients in place: self = -self.
        pub fn negate(self: *Self) void {
            negateInPlaceSimd(self.values);
        }

        /// L² inner product: ⟨self, other⟩ = Σᵢ selfᵢ · otherᵢ.
        /// This is the flat (unweighted) inner product — the Hodge-weighted
        /// version will come with the Hodge star in M2.
        pub fn inner_product(self: Self, other: Self) f64 {
            std.debug.assert(self.values.len == other.values.len);
            return innerProductSimd(self.values, other.values);
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

const Mesh2D = topology.Mesh(2, 2);

test "Cochain(Mesh, 0) has one value per vertex" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_vertices()), c.values.len);
}

test "Cochain(Mesh, 1) has one value per edge" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_edges()), c.values.len);
}

test "Cochain(Mesh, 2) has one value per face" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 2, Primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_faces()), c.values.len);
}

test "Cochain initializes to zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
    defer c.deinit(allocator);

    for (c.values) |v| {
        try testing.expectEqual(@as(f64, 0), v);
    }
}

test "Cochain degree is accessible at comptime" {
    const C0 = Cochain(Mesh2D, 0, Primal);
    const C1 = Cochain(Mesh2D, 1, Primal);
    const C2 = Cochain(Mesh2D, 2, Primal);

    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 1), C1.degree);
    try testing.expectEqual(@as(comptime_int, 2), C2.degree);
}

test "different degree cochains are distinct types" {
    // This test verifies that Cochain(Mesh, 0) and Cochain(Mesh, 1) are
    // different types — the foundation for compile-time degree enforcement.
    const C0 = Cochain(Mesh2D, 0, Primal);
    const C1 = Cochain(Mesh2D, 1, Primal);
    try testing.expect(C0 != C1);
}

test "primal and dual cochains of the same degree are distinct types" {
    const Primal1 = Cochain(Mesh2D, 1, Primal);
    const Dual1 = Cochain(Mesh2D, 1, Dual);
    try testing.expect(Primal1 != Dual1);

    // Duality is recoverable at comptime.
    try testing.expect(Primal1.duality == Primal);
    try testing.expect(Dual1.duality == Dual);
}

test "Cochain degree bounded by mesh dimension" {
    // Mesh(3, 3) supports 0-, 1-, 2-, and 3-cochains.
    const Mesh3D = topology.Mesh(3, 3);
    const C0 = Cochain(Mesh3D, 0, Primal);
    const C3 = Cochain(Mesh3D, 3, Primal);
    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 3), C3.degree);
}

// ═══════════════════════════════════════════════════════════════════════════
// Arithmetic tests
// ═══════════════════════════════════════════════════════════════════════════

test "add two 0-cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var a = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
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
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var a = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
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
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
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
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
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
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // 1×1 grid has 4 vertices
    var a = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
    defer a.deinit(allocator);
    var b = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
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

test "SIMD add kernel matches scalar reference for tail lengths" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x51ADADD0);

    const lengths = [_]usize{ 1, 2, 3, 5, 7, 9, 15, 17, 31 };
    for (lengths) |len| {
        const lhs_scalar = try allocator.alloc(f64, len);
        defer allocator.free(lhs_scalar);
        const lhs_simd = try allocator.alloc(f64, len);
        defer allocator.free(lhs_simd);
        const rhs = try allocator.alloc(f64, len);
        defer allocator.free(rhs);

        for (0..len) |i| {
            const left = rng.random().float(f64) * 200.0 - 100.0;
            const right = rng.random().float(f64) * 200.0 - 100.0;
            lhs_scalar[i] = left;
            lhs_simd[i] = left;
            rhs[i] = right;
        }

        addAssignScalar(lhs_scalar, rhs);
        addAssignSimd(lhs_simd, rhs);

        for (lhs_scalar, lhs_simd) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-14);
        }
    }
}

test "SIMD scale and negate kernels match scalar reference for tail lengths" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x51AD5CA1);

    const lengths = [_]usize{ 1, 4, 6, 8, 10, 13, 16, 19, 33 };
    for (lengths) |len| {
        const scale_scalar = try allocator.alloc(f64, len);
        defer allocator.free(scale_scalar);
        const scale_simd = try allocator.alloc(f64, len);
        defer allocator.free(scale_simd);
        const negate_scalar = try allocator.alloc(f64, len);
        defer allocator.free(negate_scalar);
        const negate_simd = try allocator.alloc(f64, len);
        defer allocator.free(negate_simd);

        const alpha = rng.random().float(f64) * 20.0 - 10.0;

        for (0..len) |i| {
            const value = rng.random().float(f64) * 200.0 - 100.0;
            scale_scalar[i] = value;
            scale_simd[i] = value;
            negate_scalar[i] = value;
            negate_simd[i] = value;
        }

        scaleInPlaceScalar(scale_scalar, alpha);
        scaleInPlaceSimd(scale_simd, alpha);
        negateInPlaceScalar(negate_scalar);
        negateInPlaceSimd(negate_simd);

        for (scale_scalar, scale_simd) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-14);
        }
        for (negate_scalar, negate_simd) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-14);
        }
    }
}

test "SIMD inner product kernel matches scalar reference for tail lengths" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x51ADD071);

    const lengths = [_]usize{ 1, 2, 5, 7, 11, 16, 18, 29, 65 };
    for (lengths) |len| {
        const lhs = try allocator.alloc(f64, len);
        defer allocator.free(lhs);
        const rhs = try allocator.alloc(f64, len);
        defer allocator.free(rhs);

        for (0..len) |i| {
            lhs[i] = rng.random().float(f64) * 200.0 - 100.0;
            rhs[i] = rng.random().float(f64) * 200.0 - 100.0;
        }

        const expected = innerProductScalar(lhs, rhs);
        const actual = innerProductSimd(lhs, rhs);
        try testing.expectApproxEqRel(expected, actual, 1e-12);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Property tests: cochain arithmetic
// ═══════════════════════════════════════════════════════════════════════════

test "add is commutative for random 0-cochains (1000 trials)" {
    // a + b = b + a for all cochains a, b.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_ADD_00);

    for (0..1000) |_| {
        var a = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer a.deinit(allocator);
        var b = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer b.deinit(allocator);

        // Save copies for the second order of addition.
        var a_copy = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer a_copy.deinit(allocator);
        var b_copy = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer b_copy.deinit(allocator);

        for (a.values, a_copy.values) |*av, *acv| {
            const val = rng.random().float(f64) * 200.0 - 100.0;
            av.* = val;
            acv.* = val;
        }
        for (b.values, b_copy.values) |*bv, *bcv| {
            const val = rng.random().float(f64) * 200.0 - 100.0;
            bv.* = val;
            bcv.* = val;
        }

        // a += b
        a.add(b);
        // b_copy += a_copy  (b + a)
        b_copy.add(a_copy);

        for (a.values, b_copy.values) |ab, ba| {
            try testing.expectApproxEqAbs(ab, ba, 1e-14);
        }
    }
}

test "add is associative for random 1-cochains (1000 trials)" {
    // (a + b) + c = a + (b + c) for all cochains a, b, c.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_A55_01);

    for (0..1000) |_| {
        var a1 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer a1.deinit(allocator);
        var a2 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer a2.deinit(allocator);
        var b1 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer b1.deinit(allocator);
        var b2 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer b2.deinit(allocator);
        var c1 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer c1.deinit(allocator);
        var c2 = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer c2.deinit(allocator);

        for (0..a1.values.len) |i| {
            const va = rng.random().float(f64) * 200.0 - 100.0;
            const vb = rng.random().float(f64) * 200.0 - 100.0;
            const vc = rng.random().float(f64) * 200.0 - 100.0;
            a1.values[i] = va;
            a2.values[i] = va;
            b1.values[i] = vb;
            b2.values[i] = vb;
            c1.values[i] = vc;
            c2.values[i] = vc;
        }

        // (a + b) + c
        a1.add(b1);
        a1.add(c1);

        // a + (b + c)
        b2.add(c2);
        a2.add(b2);

        for (a1.values, a2.values) |lhs, rhs| {
            try testing.expectApproxEqAbs(lhs, rhs, 1e-13);
        }
    }
}

test "scale distributes over add for random cochains (1000 trials)" {
    // α(a + b) = αa + αb for all cochains a, b and scalar α.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 2, 2.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_D15_02);

    for (0..1000) |_| {
        var a1 = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer a1.deinit(allocator);
        var a2 = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer a2.deinit(allocator);
        var b1 = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer b1.deinit(allocator);
        var b2 = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer b2.deinit(allocator);

        const alpha = rng.random().float(f64) * 20.0 - 10.0;

        for (0..a1.values.len) |i| {
            const va = rng.random().float(f64) * 200.0 - 100.0;
            const vb = rng.random().float(f64) * 200.0 - 100.0;
            a1.values[i] = va;
            a2.values[i] = va;
            b1.values[i] = vb;
            b2.values[i] = vb;
        }

        // α(a + b)
        a1.add(b1);
        a1.scale(alpha);

        // αa + αb
        a2.scale(alpha);
        b2.scale(alpha);
        a2.add(b2);

        for (a1.values, a2.values) |lhs, rhs| {
            try testing.expectApproxEqAbs(lhs, rhs, 1e-12);
        }
    }
}

test "negate is self-inverse for random cochains (1000 trials)" {
    // negate(negate(a)) = a for all cochains a.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_4E9_03);

    for (0..1000) |_| {
        var a = try Cochain(Mesh2D, 2, Primal).init(allocator, &mesh);
        defer a.deinit(allocator);

        const original = try allocator.alloc(f64, a.values.len);
        defer allocator.free(original);

        for (a.values, original) |*v, *o| {
            const val = rng.random().float(f64) * 200.0 - 100.0;
            v.* = val;
            o.* = val;
        }

        a.negate();
        a.negate();

        for (a.values, original) |v, o| {
            try testing.expectApproxEqAbs(o, v, 1e-15);
        }
    }
}

test "inner product is symmetric for random cochains (1000 trials)" {
    // ⟨a, b⟩ = ⟨b, a⟩ for all cochains a, b.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_1F5_04);

    for (0..1000) |_| {
        var a = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer a.deinit(allocator);
        var b = try Cochain(Mesh2D, 1, Primal).init(allocator, &mesh);
        defer b.deinit(allocator);

        for (a.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (b.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        const ab = a.inner_product(b);
        const ba = b.inner_product(a);
        try testing.expectApproxEqRel(ab, ba, 1e-14);
    }
}

test "norm_squared is non-negative for random cochains (1000 trials)" {
    // ‖a‖² ≥ 0 for all cochains a, with equality iff a = 0.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xC0C_402_05);

    for (0..1000) |_| {
        var a = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
        defer a.deinit(allocator);
        for (a.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        try testing.expect(a.norm_squared() >= 0.0);
    }

    // Zero cochain has zero norm.
    var zero = try Cochain(Mesh2D, 0, Primal).init(allocator, &mesh);
    defer zero.deinit(allocator);
    try testing.expectEqual(@as(f64, 0.0), zero.norm_squared());
}
