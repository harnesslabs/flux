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
    _ = allocator;
    _ = radius;
    _ = radial_segments;
    _ = angular_segments;
    @panic("not yet implemented");
}

pub fn torus(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    major_radius: f64,
    minor_radius: f64,
    major_segments: u32,
    minor_segments: u32,
) !MeshType {
    _ = allocator;
    _ = major_radius;
    _ = minor_radius;
    _ = major_segments;
    _ = minor_segments;
    @panic("not yet implemented");
}
