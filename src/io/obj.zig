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
    var raw_mesh = RawMesh{};
    errdefer raw_mesh.deinit(allocator);

    var line_iter = std.mem.splitScalar(u8, bytes, '\n');
    while (line_iter.next()) |line_raw| {
        const line = std.mem.trim(u8, line_raw, " \t\r");
        if (line.len == 0 or line[0] == '#') continue;

        if (std.mem.startsWith(u8, line, "v ")) {
            var tokens = std.mem.tokenizeAny(u8, line[2..], " \t");
            const x_text = tokens.next() orelse return Error.InvalidVertexRecord;
            const y_text = tokens.next() orelse return Error.InvalidVertexRecord;
            const z_text = tokens.next() orelse return Error.InvalidVertexRecord;
            const x = try std.fmt.parseFloat(f64, x_text);
            const y = try std.fmt.parseFloat(f64, y_text);
            const z = try std.fmt.parseFloat(f64, z_text);
            try raw_mesh.vertices.append(allocator, .{ x, y, z });
            continue;
        }

        if (std.mem.startsWith(u8, line, "f ")) {
            const face_start: u32 = @intCast(raw_mesh.face_vertex_indices.items.len);
            var face_len: u32 = 0;
            var tokens = std.mem.tokenizeAny(u8, line[2..], " \t");
            while (tokens.next()) |token| {
                if (token.len == 0) continue;
                var vertex_token_iter = std.mem.splitScalar(u8, token, '/');
                const vertex_token = vertex_token_iter.first();
                const vertex_index_1_based = try std.fmt.parseInt(i32, vertex_token, 10);
                if (vertex_index_1_based <= 0) {
                    if (vertex_index_1_based < 0) return Error.NegativeIndexNotSupported;
                    return Error.InvalidIndex;
                }
                const vertex_index: u32 = @intCast(vertex_index_1_based - 1);
                if (vertex_index >= raw_mesh.vertices.items.len) return Error.InvalidIndex;
                try raw_mesh.face_vertex_indices.append(allocator, vertex_index);
                face_len += 1;
            }

            if (face_len < 3) return Error.EmptyFace;
            try raw_mesh.faces.append(allocator, .{ .start = face_start, .len = face_len });
        }
    }

    return raw_mesh;
}

pub fn triangulate(allocator: std.mem.Allocator, raw_mesh: *const RawMesh) !TriangulatedSurface {
    var triangle_count: usize = 0;
    for (raw_mesh.faces.items) |face| {
        triangle_count += face.len - 2;
    }

    const triangles = try allocator.alloc([3]u32, triangle_count);
    errdefer allocator.free(triangles);

    var triangle_write: usize = 0;
    for (raw_mesh.faces.items) |face| {
        const vertex_indices = raw_mesh.face_vertex_indices.items[face.start .. face.start + face.len];
        const anchor = vertex_indices[0];
        for (1..vertex_indices.len - 1) |idx| {
            triangles[triangle_write] = .{ anchor, vertex_indices[idx], vertex_indices[idx + 1] };
            triangle_write += 1;
        }
    }
    std.debug.assert(triangle_write == triangle_count);

    return .{
        .vertex_coords = raw_mesh.vertices.items,
        .triangles = triangles,
    };
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

test "OBJ parser accepts mixed triangle and quad faces with slash indices" {
    const allocator = testing.allocator;
    const source =
        \\v 0 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 1 0
        \\v 2 0 0
        \\vt 0 0
        \\vt 1 0
        \\vt 1 1
        \\vt 0 1
        \\vt 2 0
        \\vn 0 0 1
        \\f 1/1/1 2/2/1 3/3/1
        \\f 1/1/1 3/3/1 4/4/1 5/5/1
    ;

    var raw_mesh = try parse_obj(allocator, source);
    defer raw_mesh.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), raw_mesh.faces.items.len);
    try testing.expectEqual(@as(u32, 3), raw_mesh.faces.items[0].len);
    try testing.expectEqual(@as(u32, 4), raw_mesh.faces.items[1].len);

    var triangulated = try triangulate(allocator, &raw_mesh);
    defer triangulated.deinit(allocator);

    try testing.expectEqual(@as(usize, 3), triangulated.triangles.len);
}
