const std = @import("std");

pub fn plane(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    nx: u32,
    ny: u32,
    width: f64,
    height: f64,
) !MeshType {
    return MeshType.uniform_grid(allocator, nx, ny, width, height);
}

pub fn sphere(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    radius: f64,
    refinement: u32,
) !MeshType {
    return MeshType.sphere(allocator, radius, refinement);
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

    var vertex_coords = try std.ArrayList([MeshType.embedding_dimension]f64).initCapacity(allocator, vertex_count);
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

    var vertex_coords = try std.ArrayList([MeshType.embedding_dimension]f64).initCapacity(allocator, vertex_count);
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
