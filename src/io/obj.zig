//! Wavefront `.obj` surface-mesh import.
//!
//! The parser is format-aware only: it produces a `RawMesh` preserving polygon
//! faces of any valence. A separate `triangulate()` pass converts those faces
//! into simplices so the topology layer can build a simplicial mesh.

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");

pub const Face = struct {
    start: u32,
    len: u32,
};

pub const RawMesh = struct {
    vertices: std.ArrayListUnmanaged([3]f64) = .{},
    face_vertex_indices: std.ArrayListUnmanaged(u32) = .{},
    faces: std.ArrayListUnmanaged(Face) = .{},

    pub fn deinit(self: *RawMesh, allocator: std.mem.Allocator) void {
        self.vertices.deinit(allocator);
        self.face_vertex_indices.deinit(allocator);
        self.faces.deinit(allocator);
    }
};

pub const TriangulatedSurface = struct {
    vertex_coords: []const [3]f64,
    triangles: []const [3]u32,

    pub fn deinit(self: *TriangulatedSurface, allocator: std.mem.Allocator) void {
        allocator.free(self.triangles);
    }
};

pub const Error = error{
    InvalidVertexRecord,
    InvalidFaceRecord,
    InvalidIndex,
    NegativeIndexNotSupported,
    EmptyFace,
    NotYetImplemented,
};

pub fn parse_obj(allocator: std.mem.Allocator, bytes: []const u8) !RawMesh {
    _ = allocator;
    _ = bytes;
    return Error.NotYetImplemented;
}

pub fn triangulate(allocator: std.mem.Allocator, raw_mesh: *const RawMesh) !TriangulatedSurface {
    _ = allocator;
    _ = raw_mesh;
    return Error.NotYetImplemented;
}

test "OBJ parser preserves quad face valence in RawMesh" {
    const allocator = testing.allocator;
    const source =
        \\# unit square in the z = 0 plane
        \\v 0 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 1 0
        \\f 1 2 3 4
    ;

    var raw_mesh = try parse_obj(allocator, source);
    defer raw_mesh.deinit(allocator);

    try testing.expectEqual(@as(usize, 4), raw_mesh.vertices.items.len);
    try testing.expectEqual(@as(usize, 1), raw_mesh.faces.items.len);
    try testing.expectEqual(@as(u32, 4), raw_mesh.faces.items[0].len);
}

test "OBJ triangulate splits one quad into two triangles" {
    const allocator = testing.allocator;
    const source =
        \\v 0 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 1 0
        \\f 1 2 3 4
    ;

    var raw_mesh = try parse_obj(allocator, source);
    defer raw_mesh.deinit(allocator);

    var triangulated = try triangulate(allocator, &raw_mesh);
    defer triangulated.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), triangulated.triangles.len);
}

test "OBJ import builds Mesh(3, 2) from triangulated planar quads" {
    const allocator = testing.allocator;
    const source =
        \\v 0 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 1 0
        \\v 2 0 0
        \\v 2 1 0
        \\f 1 2 3 4
        \\f 2 5 6 3
    ;

    var raw_mesh = try parse_obj(allocator, source);
    defer raw_mesh.deinit(allocator);

    var triangulated = try triangulate(allocator, &raw_mesh);
    defer triangulated.deinit(allocator);

    var mesh = try topology.Mesh(3, 2).from_triangles(allocator, triangulated.vertex_coords, triangulated.triangles);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(u32, 6), mesh.num_vertices());
    try testing.expectEqual(@as(u32, 4), mesh.num_faces());
}
