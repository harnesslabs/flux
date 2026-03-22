//! Discrete k-forms (cochains) on simplicial meshes.
//!
//! A k-cochain assigns a real value to every k-cell of a mesh. The degree `k`
//! and mesh type are comptime parameters, so the exterior derivative and other
//! operators can enforce dimensional compatibility as compile errors.

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
pub fn Cochain(comptime MeshType: type, comptime k: comptime_int) type {
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

        /// Coefficient values — one per k-cell.
        values: []f64,
        /// The mesh this cochain is defined on.
        mesh: *const MeshType,

        /// Allocate a zero-initialized cochain on the given mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            const num_cells: u32 = switch (k) {
                0 => mesh.num_vertices(),
                1 => mesh.num_edges(),
                2 => mesh.num_faces(),
                else => @compileError("unsupported degree for num_cells lookup"),
            };
            const values = try allocator.alloc(f64, num_cells);
            @memset(values, 0);
            return .{ .values = values, .mesh = mesh };
        }

        /// Free the coefficient storage.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.values);
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

    try testing.expectEqual(@as(usize, mesh.num_vertices()), c.values.len);
}

test "Cochain(Mesh, 1) has one value per edge" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 1).init(allocator, &mesh);
    defer c.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_edges()), c.values.len);
}

test "Cochain(Mesh, 2) has one value per face" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var c = try Cochain(Mesh2D, 2).init(allocator, &mesh);
    defer c.deinit(allocator);

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

    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 1), C1.degree);
    try testing.expectEqual(@as(comptime_int, 2), C2.degree);
}

test "different degree cochains are distinct types" {
    // This test verifies that Cochain(Mesh, 0) and Cochain(Mesh, 1) are
    // different types — the foundation for compile-time degree enforcement.
    const C0 = Cochain(Mesh2D, 0);
    const C1 = Cochain(Mesh2D, 1);
    try testing.expect(C0 != C1);
}

test "Cochain degree bounded by mesh dimension" {
    // Mesh(3) supports 0-, 1-, 2-, and 3-cochains.
    const Mesh3D = topology.Mesh(3);
    const C0 = Cochain(Mesh3D, 0);
    const C3 = Cochain(Mesh3D, 3);
    try testing.expectEqual(@as(comptime_int, 0), C0.degree);
    try testing.expectEqual(@as(comptime_int, 3), C3.degree);
}
