const std = @import("std");

pub fn plane(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    nx: u32,
    ny: u32,
    width: f64,
    height: f64,
) !MeshType {
    comptime {
        if (MeshType.embedding_dimension != 2 or MeshType.topological_dimension != 2) {
            @compileError("plane is only available for Mesh(2, 2)");
        }
    }
    std.debug.assert(nx > 0);
    std.debug.assert(ny > 0);
    std.debug.assert(width > 0.0);
    std.debug.assert(height > 0.0);

    const vertex_count: usize = (nx + 1) * (ny + 1);
    const face_count: usize = 2 * nx * ny;
    var vertex_coords = try std.ArrayList([2]f64).initCapacity(allocator, vertex_count);
    defer vertex_coords.deinit(allocator);
    var face_vertices = try std.ArrayList([3]u32).initCapacity(allocator, face_count);
    defer face_vertices.deinit(allocator);

    const dx = width / @as(f64, @floatFromInt(nx));
    const dy = height / @as(f64, @floatFromInt(ny));

    for (0..nx + 1) |i_usize| {
        const i: u32 = @intCast(i_usize);
        for (0..ny + 1) |j_usize| {
            const j: u32 = @intCast(j_usize);
            vertex_coords.appendAssumeCapacity(.{
                @as(f64, @floatFromInt(i)) * dx,
                @as(f64, @floatFromInt(j)) * dy,
            });
        }
    }

    for (0..nx) |i_usize| {
        const i: u32 = @intCast(i_usize);
        for (0..ny) |j_usize| {
            const j: u32 = @intCast(j_usize);
            const sw = planeVertexIndex(i, j, ny);
            const se = planeVertexIndex(i + 1, j, ny);
            const nw = planeVertexIndex(i, j + 1, ny);
            const ne = planeVertexIndex(i + 1, j + 1, ny);
            face_vertices.appendAssumeCapacity(.{ sw, se, ne });
            face_vertices.appendAssumeCapacity(.{ sw, ne, nw });
        }
    }

    return MeshType.from_triangles(allocator, vertex_coords.items, face_vertices.items);
}

pub fn sphere(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    radius: f64,
    refinement: u32,
) !MeshType {
    comptime {
        if (MeshType.embedding_dimension != 3 or MeshType.topological_dimension != 2) {
            @compileError("sphere is only available for Mesh(3, 2)");
        }
    }
    std.debug.assert(radius > 0.0);

    var polyhedron = try buildRefinedSphereOctahedron(allocator, radius, refinement);
    defer polyhedron.deinit(allocator);

    const oriented_faces = try allocator.dupe([3]u32, polyhedron.faces);
    defer allocator.free(oriented_faces);
    orientSphereFacesOutward(polyhedron.vertices, oriented_faces);

    return MeshType.from_triangles(allocator, polyhedron.vertices, oriented_faces);
}

pub fn disk(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    radius: f64,
    radial_segments: u32,
    angular_segments: u32,
) !MeshType {
    comptime {
        if (MeshType.embedding_dimension != 2 or MeshType.topological_dimension != 2) {
            @compileError("disk is only available for Mesh(2, 2)");
        }
    }
    std.debug.assert(radius > 0.0);
    std.debug.assert(radial_segments > 0);
    std.debug.assert(angular_segments >= 3);

    const vertex_count: usize = 1 + radial_segments * angular_segments;
    const face_count: usize = angular_segments + (radial_segments - 1) * angular_segments * 2;

    var vertex_coords = try std.ArrayList([2]f64).initCapacity(allocator, vertex_count);
    defer vertex_coords.deinit(allocator);

    var face_vertices = try std.ArrayList([3]u32).initCapacity(allocator, face_count);
    defer face_vertices.deinit(allocator);

    vertex_coords.appendAssumeCapacity(.{ 0.0, 0.0 });
    for (1..radial_segments + 1) |ring_usize| {
        const ring: u32 = @intCast(ring_usize);
        const ring_radius = radius * @as(f64, @floatFromInt(ring)) / @as(f64, @floatFromInt(radial_segments));
        for (0..angular_segments) |sector_usize| {
            const sector: u32 = @intCast(sector_usize);
            const theta = 2.0 * std.math.pi * @as(f64, @floatFromInt(sector)) / @as(f64, @floatFromInt(angular_segments));
            vertex_coords.appendAssumeCapacity(.{
                ring_radius * @cos(theta),
                ring_radius * @sin(theta),
            });
        }
    }

    for (0..angular_segments) |sector_usize| {
        const sector: u32 = @intCast(sector_usize);
        face_vertices.appendAssumeCapacity(.{
            0,
            diskVertexIndex(1, sector, angular_segments),
            diskVertexIndex(1, (sector + 1) % angular_segments, angular_segments),
        });
    }

    for (1..radial_segments) |ring_usize| {
        const ring: u32 = @intCast(ring_usize);
        for (0..angular_segments) |sector_usize| {
            const sector: u32 = @intCast(sector_usize);
            const inner_0 = diskVertexIndex(ring, sector, angular_segments);
            const inner_1 = diskVertexIndex(ring, (sector + 1) % angular_segments, angular_segments);
            const outer_0 = diskVertexIndex(ring + 1, sector, angular_segments);
            const outer_1 = diskVertexIndex(ring + 1, (sector + 1) % angular_segments, angular_segments);

            face_vertices.appendAssumeCapacity(.{ inner_0, outer_0, outer_1 });
            face_vertices.appendAssumeCapacity(.{ inner_0, outer_1, inner_1 });
        }
    }

    return MeshType.from_triangles(allocator, vertex_coords.items, face_vertices.items);
}

pub fn torus(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    major_radius: f64,
    minor_radius: f64,
    major_segments: u32,
    minor_segments: u32,
) !MeshType {
    comptime {
        if (MeshType.embedding_dimension != 3 or MeshType.topological_dimension != 2) {
            @compileError("torus is only available for Mesh(3, 2)");
        }
    }
    std.debug.assert(major_radius > 0.0);
    std.debug.assert(minor_radius > 0.0);
    std.debug.assert(major_radius > minor_radius);
    std.debug.assert(major_segments >= 3);
    std.debug.assert(minor_segments >= 3);

    const vertex_count: usize = major_segments * minor_segments;
    const face_count: usize = 2 * major_segments * minor_segments;

    var vertex_coords = try std.ArrayList([3]f64).initCapacity(allocator, vertex_count);
    defer vertex_coords.deinit(allocator);

    var face_vertices = try std.ArrayList([3]u32).initCapacity(allocator, face_count);
    defer face_vertices.deinit(allocator);

    for (0..major_segments) |major_usize| {
        const major: u32 = @intCast(major_usize);
        const u = 2.0 * std.math.pi * @as(f64, @floatFromInt(major)) / @as(f64, @floatFromInt(major_segments));
        const cos_u = @cos(u);
        const sin_u = @sin(u);

        for (0..minor_segments) |minor_usize| {
            const minor: u32 = @intCast(minor_usize);
            const v = 2.0 * std.math.pi * @as(f64, @floatFromInt(minor)) / @as(f64, @floatFromInt(minor_segments));
            const cos_v = @cos(v);
            const sin_v = @sin(v);
            const ring_radius = major_radius + minor_radius * cos_v;
            vertex_coords.appendAssumeCapacity(.{
                ring_radius * cos_u,
                ring_radius * sin_u,
                minor_radius * sin_v,
            });
        }
    }

    for (0..major_segments) |major_usize| {
        const major: u32 = @intCast(major_usize);
        const major_next = (major + 1) % major_segments;
        for (0..minor_segments) |minor_usize| {
            const minor: u32 = @intCast(minor_usize);
            const minor_next = (minor + 1) % minor_segments;
            const v00 = torusVertexIndex(major, minor, minor_segments);
            const v10 = torusVertexIndex(major_next, minor, minor_segments);
            const v11 = torusVertexIndex(major_next, minor_next, minor_segments);
            const v01 = torusVertexIndex(major, minor_next, minor_segments);

            face_vertices.appendAssumeCapacity(.{ v00, v10, v11 });
            face_vertices.appendAssumeCapacity(.{ v00, v11, v01 });
        }
    }

    return MeshType.from_triangles(allocator, vertex_coords.items, face_vertices.items);
}

fn planeVertexIndex(i: u32, j: u32, ny: u32) u32 {
    return i * (ny + 1) + j;
}

fn diskVertexIndex(ring: u32, sector: u32, angular_segments: u32) u32 {
    std.debug.assert(ring >= 1);
    return 1 + (ring - 1) * angular_segments + sector;
}

fn torusVertexIndex(major: u32, minor: u32, minor_segments: u32) u32 {
    return major * minor_segments + minor;
}

const SpherePolyhedron = struct {
    vertices: [][3]f64,
    faces: [][3]u32,

    fn deinit(self: *SpherePolyhedron, allocator: std.mem.Allocator) void {
        allocator.free(self.faces);
        allocator.free(self.vertices);
    }
};

fn buildRefinedSphereOctahedron(allocator: std.mem.Allocator, radius: f64, refinement: u32) !SpherePolyhedron {
    var polyhedron = SpherePolyhedron{
        .vertices = try allocator.dupe([3]f64, &initialSphereOctahedronVertices),
        .faces = try allocator.dupe([3]u32, &initialSphereOctahedronFaces),
    };
    errdefer polyhedron.deinit(allocator);

    scaleSphereVertices(polyhedron.vertices, radius);
    for (0..refinement) |_| {
        const next = try refineSpherePolyhedron(allocator, polyhedron, radius);
        polyhedron.deinit(allocator);
        polyhedron = next;
    }
    return polyhedron;
}

fn refineSpherePolyhedron(allocator: std.mem.Allocator, polyhedron: SpherePolyhedron, radius: f64) !SpherePolyhedron {
    var vertices = try std.ArrayList([3]f64).initCapacity(allocator, polyhedron.vertices.len + polyhedron.faces.len);
    defer vertices.deinit(allocator);
    try vertices.appendSlice(allocator, polyhedron.vertices);

    var faces = try std.ArrayList([3]u32).initCapacity(allocator, polyhedron.faces.len * 4);
    defer faces.deinit(allocator);

    var midpoint_map = std.AutoHashMap([2]u32, u32).init(allocator);
    defer midpoint_map.deinit();

    for (polyhedron.faces) |face| {
        const ab = try sphereMidpointIndex(allocator, &vertices, &midpoint_map, face[0], face[1], radius);
        const bc = try sphereMidpointIndex(allocator, &vertices, &midpoint_map, face[1], face[2], radius);
        const ca = try sphereMidpointIndex(allocator, &vertices, &midpoint_map, face[2], face[0], radius);

        try faces.append(allocator, .{ face[0], ab, ca });
        try faces.append(allocator, .{ face[1], bc, ab });
        try faces.append(allocator, .{ face[2], ca, bc });
        try faces.append(allocator, .{ ab, bc, ca });
    }

    return .{
        .vertices = try vertices.toOwnedSlice(allocator),
        .faces = try faces.toOwnedSlice(allocator),
    };
}

fn sphereMidpointIndex(
    allocator: std.mem.Allocator,
    vertices: *std.ArrayList([3]f64),
    midpoint_map: *std.AutoHashMap([2]u32, u32),
    a: u32,
    b: u32,
    radius: f64,
) !u32 {
    const key = canonicalEdgeKey(a, b);
    const gop = try midpoint_map.getOrPut(key);
    if (gop.found_existing) return gop.value_ptr.*;

    const midpoint = normalizeSpherePoint(scaleSpherePoint(addSpherePoints(vertices.items[a], vertices.items[b]), 0.5), radius);
    const new_index: u32 = @intCast(vertices.items.len);
    try vertices.append(allocator, midpoint);
    gop.value_ptr.* = new_index;
    return new_index;
}

fn orientSphereFacesOutward(embedded_vertices: []const [3]f64, faces: [][3]u32) void {
    for (faces) |*face| {
        const a = embedded_vertices[face.*[0]];
        const b = embedded_vertices[face.*[1]];
        const c = embedded_vertices[face.*[2]];
        const normal = cross3(subSpherePoints(b, a), subSpherePoints(c, a));
        const centroid = normalizeSpherePoint(scaleSpherePoint(addSpherePoints(addSpherePoints(a, b), c), 1.0 / 3.0), 1.0);
        if (dot3(normal, centroid) < 0.0) {
            std.mem.swap(u32, &face.*[1], &face.*[2]);
        }
    }
}

fn scaleSphereVertices(vertices: [][3]f64, radius: f64) void {
    for (vertices) |*vertex| {
        vertex.* = scaleSpherePoint(vertex.*, radius);
    }
}

fn addSpherePoints(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

fn subSpherePoints(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

fn scaleSpherePoint(a: [3]f64, scale: f64) [3]f64 {
    return .{ a[0] * scale, a[1] * scale, a[2] * scale };
}

fn dot3(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

fn norm3(a: [3]f64) f64 {
    return @sqrt(dot3(a, a));
}

fn normalizeSpherePoint(a: [3]f64, radius: f64) [3]f64 {
    return scaleSpherePoint(a, radius / norm3(a));
}

fn cross3(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

fn canonicalEdgeKey(a: u32, b: u32) [2]u32 {
    return if (a < b) .{ a, b } else .{ b, a };
}

const initialSphereOctahedronVertices = [_][3]f64{
    .{ 1.0, 0.0, 0.0 },
    .{ -1.0, 0.0, 0.0 },
    .{ 0.0, 1.0, 0.0 },
    .{ 0.0, -1.0, 0.0 },
    .{ 0.0, 0.0, 1.0 },
    .{ 0.0, 0.0, -1.0 },
};

const initialSphereOctahedronFaces = [_][3]u32{
    .{ 0, 2, 4 },
    .{ 2, 1, 4 },
    .{ 1, 3, 4 },
    .{ 3, 0, 4 },
    .{ 2, 0, 5 },
    .{ 1, 2, 5 },
    .{ 3, 1, 5 },
    .{ 0, 3, 5 },
};
