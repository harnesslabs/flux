const std = @import("std");
const two_dimensional = @import("two_dimensional.zig");
const three_dimensional = @import("three_dimensional.zig");

pub fn Mesh(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => two_dimensional.Mesh,
        3 => three_dimensional.Mesh,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn Config(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => two_dimensional.ConfigImpl,
        3 => three_dimensional.ConfigImpl,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn RunResult(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => two_dimensional.RunResultImpl,
        3 => three_dimensional.RunResultImpl,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn State(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => two_dimensional.StateImpl,
        3 => three_dimensional.StateImpl,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn run(comptime topological_dimension: u8, allocator: std.mem.Allocator, config: Config(topological_dimension), writer: anytype) !RunResult(topological_dimension) {
    return switch (topological_dimension) {
        2 => try two_dimensional.runImpl(allocator, config, writer),
        3 => try three_dimensional.runImpl(allocator, config, writer),
        else => unreachable,
    };
}

pub fn step(comptime topological_dimension: u8, allocator: std.mem.Allocator, state: *State(topological_dimension), dt: f64) !void {
    switch (topological_dimension) {
        2 => try two_dimensional.stepImpl(allocator, state, dt),
        3 => try three_dimensional.stepImpl(allocator, state, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(comptime topological_dimension: u8, allocator: std.mem.Allocator, state: *State(topological_dimension)) !void {
    switch (topological_dimension) {
        2 => try two_dimensional.seedReferenceMode(allocator, state),
        3 => try three_dimensional.seedReferenceMode(allocator, state),
        else => unreachable,
    }
}

pub fn conservedQuantity(comptime topological_dimension: u8, allocator: std.mem.Allocator, state: *const State(topological_dimension)) !f64 {
    return switch (topological_dimension) {
        2 => two_dimensional.conservedQuantity(state),
        3 => try three_dimensional.conservedQuantity(allocator, state),
        else => unreachable,
    };
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
