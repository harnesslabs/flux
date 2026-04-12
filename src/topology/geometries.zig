const std = @import("std");
const sparse = @import("../math/sparse.zig");

pub fn plane(
    comptime MeshType: type,
    comptime Api: type,
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

    const dx = width / @as(f64, @floatFromInt(nx));
    const dy = height / @as(f64, @floatFromInt(ny));

    const vertex_count: u32 = (nx + 1) * (ny + 1);
    const horizontal_edge_count: u32 = nx * (ny + 1);
    const vertical_edge_count: u32 = (nx + 1) * ny;
    const diagonal_edge_count: u32 = nx * ny;
    const edge_count: u32 = horizontal_edge_count + vertical_edge_count + diagonal_edge_count;
    const face_count: u32 = 2 * nx * ny;

    var vertices: Api.VerticesStorage = .{};
    try vertices.ensureTotalCapacity(allocator, vertex_count);
    errdefer vertices.deinit(allocator);

    var edges_list: Api.EdgeStorage = .{};
    try edges_list.ensureTotalCapacity(allocator, edge_count);
    errdefer edges_list.deinit(allocator);

    var faces_list: Api.FaceStorage = .{};
    try faces_list.ensureTotalCapacity(allocator, face_count);
    errdefer faces_list.deinit(allocator);

    for (0..nx + 1) |i_u| {
        const fi: f64 = @floatFromInt(i_u);
        for (0..ny + 1) |j_u| {
            const fj: f64 = @floatFromInt(j_u);
            var coords: [MeshType.embedding_dimension]f64 = @splat(0);
            coords[0] = fi * dx;
            coords[1] = fj * dy;
            vertices.appendAssumeCapacity(.{ .coords = coords, .dual_volume = 0 });
        }
    }

    for (0..ny + 1) |j_u| {
        for (0..nx) |i_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            edges_list.appendAssumeCapacity(.{
                .vertices = .{ Api.vertexIndex(i, j, ny), Api.vertexIndex(i + 1, j, ny) },
                .volume = 0,
                .barycenter = @splat(0),
            });
        }
    }
    for (0..nx + 1) |i_u| {
        for (0..ny) |j_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            edges_list.appendAssumeCapacity(.{
                .vertices = .{ Api.vertexIndex(i, j, ny), Api.vertexIndex(i, j + 1, ny) },
                .volume = 0,
                .barycenter = @splat(0),
            });
        }
    }
    for (0..nx) |i_u| {
        for (0..ny) |j_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            edges_list.appendAssumeCapacity(.{
                .vertices = .{ Api.vertexIndex(i, j, ny), Api.vertexIndex(i + 1, j + 1, ny) },
                .volume = 0,
                .barycenter = @splat(0),
            });
        }
    }

    for (0..nx) |i_u| {
        const i: u32 = @intCast(i_u);
        for (0..ny) |j_u| {
            const j: u32 = @intCast(j_u);
            const sw = Api.vertexIndex(i, j, ny);
            const se = Api.vertexIndex(i + 1, j, ny);
            const nw = Api.vertexIndex(i, j + 1, ny);
            const ne = Api.vertexIndex(i + 1, j + 1, ny);
            const zero_barycenter: [MeshType.embedding_dimension]f64 = @splat(0);

            faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, se, ne }, .volume = 0, .barycenter = zero_barycenter });
            faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, ne, nw }, .volume = 0, .barycenter = zero_barycenter });
        }
    }

    var boundary_1: sparse.PackedIncidenceMatrix = undefined;
    {
        var dense = try sparse.CsrMatrix(i8).init(allocator, edge_count, vertex_count, 2 * edge_count);
        errdefer dense.deinit(allocator);
        const edge_verts = edges_list.slice().items(.vertices);
        for (0..edge_count) |edge_idx| {
            dense.row_ptr[edge_idx] = @intCast(2 * edge_idx);
            dense.col_idx[2 * edge_idx] = edge_verts[edge_idx][0];
            dense.values[2 * edge_idx] = -1;
            dense.col_idx[2 * edge_idx + 1] = edge_verts[edge_idx][1];
            dense.values[2 * edge_idx + 1] = 1;
        }
        dense.row_ptr[edge_count] = 2 * edge_count;

        boundary_1 = try sparse.PackedIncidenceMatrix.fromBoundaryCsr(allocator, 1, dense);
        dense.deinit(allocator);
    }
    errdefer boundary_1.deinit(allocator);

    var boundary_2: sparse.PackedIncidenceMatrix = undefined;
    {
        var dense = try sparse.CsrMatrix(i8).init(allocator, face_count, edge_count, 3 * face_count);
        errdefer dense.deinit(allocator);
        var face_idx: u32 = 0;
        for (0..nx) |i_u| {
            const i: u32 = @intCast(i_u);
            for (0..ny) |j_u| {
                const j: u32 = @intCast(j_u);
                const h_ij = Api.horizontalEdgeIndex(i, j, nx);
                const h_i_jp1 = Api.horizontalEdgeIndex(i, j + 1, nx);
                const v_ip1_j = Api.verticalEdgeIndex(i + 1, j, ny, horizontal_edge_count);
                const v_i_j = Api.verticalEdgeIndex(i, j, ny, horizontal_edge_count);
                const d_ij = Api.diagonalEdgeIndex(i, j, ny, horizontal_edge_count, vertical_edge_count);

                dense.row_ptr[face_idx] = 3 * face_idx;
                dense.col_idx[3 * face_idx] = h_ij;
                dense.values[3 * face_idx] = 1;
                dense.col_idx[3 * face_idx + 1] = v_ip1_j;
                dense.values[3 * face_idx + 1] = 1;
                dense.col_idx[3 * face_idx + 2] = d_ij;
                dense.values[3 * face_idx + 2] = -1;
                face_idx += 1;

                dense.row_ptr[face_idx] = 3 * face_idx;
                dense.col_idx[3 * face_idx] = h_i_jp1;
                dense.values[3 * face_idx] = -1;
                dense.col_idx[3 * face_idx + 1] = v_i_j;
                dense.values[3 * face_idx + 1] = -1;
                dense.col_idx[3 * face_idx + 2] = d_ij;
                dense.values[3 * face_idx + 2] = 1;
                face_idx += 1;
            }
        }
        dense.row_ptr[face_count] = 3 * face_count;

        boundary_2 = try sparse.PackedIncidenceMatrix.fromBoundaryCsr(allocator, 2, dense);
        dense.deinit(allocator);
    }
    errdefer boundary_2.deinit(allocator);

    const coords = vertices.slice().items(.coords);

    {
        const edge_slice = edges_list.slice();
        const edge_verts = edge_slice.items(.vertices);
        const edge_volumes = edge_slice.items(.volume);
        const edge_barycenters = edge_slice.items(.barycenter);
        for (0..edge_count) |edge_idx| {
            edge_volumes[edge_idx] = Api.euclideanDistance(coords[edge_verts[edge_idx][0]], coords[edge_verts[edge_idx][1]]);
            edge_barycenters[edge_idx] = Api.pointMidpoint(coords[edge_verts[edge_idx][0]], coords[edge_verts[edge_idx][1]]);
        }
    }

    {
        const face_slice = faces_list.slice();
        const face_verts = face_slice.items(.vertices);
        const face_volumes = face_slice.items(.volume);
        const face_barycenters = face_slice.items(.barycenter);
        for (0..face_count) |face_idx| {
            const p0 = coords[face_verts[face_idx][0]];
            const p1 = coords[face_verts[face_idx][1]];
            const p2 = coords[face_verts[face_idx][2]];
            face_volumes[face_idx] = Api.triangleArea(p0, p1, p2);
            face_barycenters[face_idx] = Api.triangleBarycenter(p0, p1, p2);
        }
    }

    const dual_edge_volumes = try allocator.alloc(f64, edge_count);
    errdefer allocator.free(dual_edge_volumes);
    var boundary_edges: []u32 = &.{};
    {
        const scratch = try allocator.alloc(u32, 3 * edge_count);
        defer allocator.free(scratch);
        const edge_face_count = scratch[0..edge_count];
        const edge_face_0 = scratch[edge_count .. 2 * edge_count];
        const edge_face_1 = scratch[2 * edge_count .. 3 * edge_count];
        @memset(edge_face_count, 0);

        for (0..face_count) |face_idx| {
            const face_edges = boundary_2.row(@intCast(face_idx));
            for (face_edges.cols) |edge_idx| {
                const count = edge_face_count[edge_idx];
                if (count == 0) {
                    edge_face_0[edge_idx] = @intCast(face_idx);
                } else {
                    edge_face_1[edge_idx] = @intCast(face_idx);
                }
                edge_face_count[edge_idx] = count + 1;
            }
        }

        const barycenters = faces_list.slice().items(.barycenter);
        const edge_verts = edges_list.slice().items(.vertices);
        var boundary_count: u32 = 0;
        for (0..edge_count) |edge_idx| {
            if (edge_face_count[edge_idx] == 2) {
                dual_edge_volumes[edge_idx] = Api.euclideanDistance(barycenters[edge_face_0[edge_idx]], barycenters[edge_face_1[edge_idx]]);
            } else if (edge_face_count[edge_idx] == 1) {
                const midpoint = Api.pointMidpoint(coords[edge_verts[edge_idx][0]], coords[edge_verts[edge_idx][1]]);
                dual_edge_volumes[edge_idx] = Api.euclideanDistance(barycenters[edge_face_0[edge_idx]], midpoint);
                boundary_count += 1;
            } else {
                return error.NonManifoldEdge;
            }
        }

        boundary_edges = try allocator.alloc(u32, boundary_count);
        var boundary_write: u32 = 0;
        for (0..edge_count) |edge_idx| {
            if (edge_face_count[edge_idx] != 1) continue;
            boundary_edges[boundary_write] = @intCast(edge_idx);
            boundary_write += 1;
        }
        std.debug.assert(boundary_write == boundary_count);
    }

    {
        const dual_volumes = vertices.slice().items(.dual_volume);
        @memset(dual_volumes, 0);

        const face_verts = faces_list.slice().items(.vertices);
        for (0..face_count) |face_idx| {
            const vs = face_verts[face_idx];
            const p0 = coords[vs[0]];
            const p1 = coords[vs[1]];
            const p2 = coords[vs[2]];

            const l01_sq = Api.distanceSquared(p0, p1);
            const l12_sq = Api.distanceSquared(p1, p2);
            const l20_sq = Api.distanceSquared(p2, p0);

            const area_4 = 4.0 * Api.triangleArea(p0, p1, p2);
            if (!(area_4 > 0.0)) return error.DegenerateTriangle;

            const cot0 = (l01_sq + l20_sq - l12_sq) / area_4;
            const cot1 = (l01_sq + l12_sq - l20_sq) / area_4;
            const cot2 = (l12_sq + l20_sq - l01_sq) / area_4;

            dual_volumes[vs[0]] += (cot1 * l20_sq + cot2 * l01_sq) / 8.0;
            dual_volumes[vs[1]] += (cot2 * l01_sq + cot0 * l12_sq) / 8.0;
            dual_volumes[vs[2]] += (cot0 * l12_sq + cot1 * l20_sq) / 8.0;
        }
    }

    var mesh = MeshType{
        .vertices = vertices,
        .simplex_lists = .{ edges_list, faces_list },
        .boundaries = .{ boundary_1, boundary_2 },
        .dual_edge_volumes = dual_edge_volumes,
        .boundary_edges = boundary_edges,
        .whitney_operators = undefined,
    };
    errdefer {
        allocator.free(mesh.boundary_edges);
        allocator.free(mesh.dual_edge_volumes);
        mesh.vertices.deinit(allocator);
        inline for (0..MeshType.topological_dimension) |simplex_idx| {
            mesh.simplex_lists[simplex_idx].deinit(allocator);
        }
        for (&mesh.boundaries) |*boundary| {
            boundary.deinit(allocator);
        }
    }

    mesh.whitney_operators = try Api.assembleWhitneyOperators(allocator, &mesh);
    return mesh;
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
