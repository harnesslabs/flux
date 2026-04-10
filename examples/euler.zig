const std = @import("std");

const euler_2d = @import("euler_2d/euler.zig");
const euler_3d = @import("euler_3d/euler.zig");

pub fn Mesh(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.Mesh2D,
        3 => euler_3d.Mesh3D,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn Config(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.Config,
        3 => euler_3d.Config,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn RunResult(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.RunResult,
        3 => euler_3d.RunResult,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn State(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.State,
        3 => euler_3d.State,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn run(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    config: Config(topological_dimension),
    writer: anytype,
) !RunResult(topological_dimension) {
    return switch (topological_dimension) {
        2 => try euler_2d.run(allocator, config, writer),
        3 => try euler_3d.run(allocator, config, writer),
        else => unreachable,
    };
}

pub fn step(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    state: *State(topological_dimension),
    dt: f64,
) !void {
    switch (topological_dimension) {
        2 => try euler_2d.step(allocator, state, dt),
        3 => try euler_3d.step(allocator, state, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    state: *State(topological_dimension),
) !void {
    if (topological_dimension == 2) {
        euler_2d.initializeVortexDipole(state);
        return;
    }
    if (topological_dimension == 3) {
        return euler_3d.seedReferenceMode(allocator, state);
    }
    unreachable;
}

pub fn conservedQuantity(
    comptime topological_dimension: u8,
    allocator: std.mem.Allocator,
    state: *const State(topological_dimension),
) !f64 {
    if (topological_dimension == 2) {
        return euler_2d.totalCirculation(state);
    }
    if (topological_dimension == 3) {
        return euler_3d.computeHelicity(allocator, state);
    }
    unreachable;
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
