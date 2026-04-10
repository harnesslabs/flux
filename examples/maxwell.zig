const std = @import("std");

const maxwell_2d = @import("maxwell_2d/maxwell.zig");
const maxwell_3d = @import("maxwell_3d/maxwell.zig");

pub fn Mesh(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => maxwell_2d.Mesh(2),
        3 => maxwell_3d.Mesh3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn Config(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => maxwell_2d.Config,
        3 => maxwell_3d.Config,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn State(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => maxwell_2d.State(Mesh(2)),
        3 => maxwell_3d.State(Mesh(3)),
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn runDriver(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    config: Config(topological_dimension),
) !void {
    switch (topological_dimension) {
        2 => try maxwell_2d.runDriver(allocator, config),
        3 => try maxwell_3d.runDriver(allocator, config),
        else => unreachable,
    }
}

pub fn makeMesh(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    config: Config(topological_dimension),
) !Mesh(topological_dimension) {
    return switch (topological_dimension) {
        2 => try Mesh(2).uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain),
        3 => try maxwell_3d.makeCavityMesh(allocator, config),
        else => unreachable,
    };
}

pub fn leapfrogStep(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    state: *State(topological_dimension),
    dt: f64,
) !void {
    switch (topological_dimension) {
        2 => try maxwell_2d.leapfrog_step(allocator, state, dt),
        3 => try maxwell_3d.leapfrogStep(allocator, state, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    state: *State(topological_dimension),
    dt: f64,
    width: f64,
    height: f64,
) !void {
    switch (topological_dimension) {
        2 => maxwell_2d.project_te10_b(state.mesh, state.B.values, -dt / 2.0, width),
        3 => try maxwell_3d.seedTm110Mode(allocator, state, dt, width, height),
        else => unreachable,
    }
}

pub fn divergenceNorm(
    allocator: std.mem.Allocator,
    state: *const State(3),
) !f64 {
    return maxwell_3d.divergenceNorm(allocator, state);
}

pub fn projectReferenceMagneticField(
    comptime topological_dimension: u8,
    mesh: *const Mesh(topological_dimension),
    values: []f64,
    t: f64,
    width: f64,
    height: f64,
) void {
    switch (topological_dimension) {
        2 => maxwell_2d.project_te10_b(mesh, values, t, width),
        3 => maxwell_3d.project_tm110_b(mesh, values, t, width, height),
        else => unreachable,
    }
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
