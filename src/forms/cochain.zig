const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");

/// A discrete k-form (cochain) on a simplicial mesh.
///
/// In DEC, a k-cochain assigns a real value to every k-cell of the mesh:
///   - 0-cochain: one value per vertex   (discrete scalar field)
///   - 1-cochain: one value per edge     (discrete circulation/flux)
///   - 2-cochain: one value per face     (discrete flux/density)
///
/// The degree `k` and mesh type are fixed at comptime, so that operators
/// like the exterior derivative `d` can enforce degree compatibility as
/// compile errors rather than runtime checks.
pub fn Cochain(comptime MeshType: type, comptime k: u2) type {
    // Validate that MeshType is actually a Mesh by checking for the
    // dimension constant and required entity accessors.
    comptime {
        if (!@hasDecl(MeshType, "dimension")) {
            @compileError("Cochain requires a Mesh type with a 'dimension' declaration");
        }
        if (k > 2) {
            @compileError("Cochain degree must be 0, 1, or 2 for a 2D mesh");
        }
    }

    return struct {
        const Self = @This();

        /// The mesh type this cochain is defined on.
        pub const MeshT = MeshType;
        /// The degree of this cochain (0 = vertices, 1 = edges, 2 = faces).
        pub const degree = k;

        /// Coefficient values — one per k-cell.
        values: []f64,
        /// The mesh this cochain is defined on (pointer to avoid copies).
        mesh: *const MeshType,

        /// Allocate a zero-initialized cochain on the given mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            const count = num_cells(mesh);
            const values = try allocator.alloc(f64, count);
            @memset(values, 0);
            return .{ .values = values, .mesh = mesh };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.values);
        }

        /// Number of k-cells in the mesh (and thus the length of this cochain).
        pub fn len(self: Self) u32 {
            return num_cells(self.mesh);
        }

        /// Return the number of k-cells for the given mesh.
        fn num_cells(mesh: *const MeshType) u32 {
            return switch (k) {
                0 => mesh.num_vertices(),
                1 => mesh.num_edges(),
                2 => mesh.num_faces(),
                else => unreachable,
            };
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

    var c = try Cochain(Mesh2D, 0).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(mesh.num_vertices(), c.len());
    try testing.expectEqual(@as(usize, mesh.num_vertices()), c.values.len);
}

test "Cochain(Mesh, 1) has one value per edge" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(mesh.num_edges(), c.len());
    try testing.expectEqual(@as(usize, mesh.num_edges()), c.values.len);
}

test "Cochain(Mesh, 2) has one value per face" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 2).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(mesh.num_faces(), c.len());
    try testing.expectEqual(@as(usize, mesh.num_faces()), c.values.len);
}

test "Cochain initializes to zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1).init(allocator, &mesh);
    defer c.deinit(allocator);

    for (c.values) |v| {
        try testing.expectEqual(@as(f64, 0), v);
    }
}

test "Cochain degree is accessible at comptime" {
    const C0 = Cochain(Mesh2D, 0);
    const C1 = Cochain(Mesh2D, 1);
    const C2 = Cochain(Mesh2D, 2);

    try testing.expectEqual(@as(u2, 0), C0.degree);
    try testing.expectEqual(@as(u2, 1), C1.degree);
    try testing.expectEqual(@as(u2, 2), C2.degree);
}

test "different degree cochains are distinct types" {
    // This test verifies that Cochain(Mesh, 0) and Cochain(Mesh, 1) are
    // different types — the foundation for compile-time degree enforcement.
    const C0 = Cochain(Mesh2D, 0);
    const C1 = Cochain(Mesh2D, 1);
    try testing.expect(C0 != C1);
}

test "Cochain works with Mesh(3) embedding" {
    const Mesh3D = topology.Mesh(3);
    const C0 = Cochain(Mesh3D, 0);
    try testing.expectEqual(@as(u2, 0), C0.degree);
}
