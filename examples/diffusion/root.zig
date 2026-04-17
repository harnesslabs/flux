const std = @import("std");
const plane = @import("plane.zig");
const sphere = @import("sphere.zig");

pub const SurfaceKind = enum {
    plane,
    sphere,
};

pub fn Mesh(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => plane.Mesh2D,
        .sphere => sphere.SurfaceMesh,
    };
}

pub fn System(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => plane.SystemImpl,
        .sphere => sphere.SystemImpl,
    };
}

pub fn Config(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => plane.ConfigImpl,
        .sphere => sphere.ConfigImpl,
    };
}

pub fn RunResult(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => plane.RunResultImpl,
        .sphere => sphere.RunResultImpl,
    };
}

pub fn ConvergenceResult(comptime surface_kind: SurfaceKind) type {
    return switch (surface_kind) {
        .plane => plane.ConvergenceResultImpl,
        .sphere => sphere.ConvergenceResultImpl,
    };
}

pub fn run(comptime surface_kind: SurfaceKind, allocator: std.mem.Allocator, config: Config(surface_kind), writer: anytype) !RunResult(surface_kind) {
    return switch (surface_kind) {
        .plane => try plane.runImpl(allocator, config, writer),
        .sphere => try sphere.runImpl(allocator, config, writer),
    };
}

pub fn runConvergenceStudy(comptime surface_kind: SurfaceKind, allocator: std.mem.Allocator, params: []const u32) ![]ConvergenceResult(surface_kind) {
    return switch (surface_kind) {
        .plane => try plane.runConvergenceStudyImpl(allocator, params),
        .sphere => try sphere.runConvergenceStudyImpl(allocator, params),
    };
}

test {
    _ = @import("tests.zig");
}
