//! Generic leapfrog (Störmer-Verlet) method family.
//!
//! The method does not own runtime state. Instead, the chosen system exposes a
//! leapfrog-compatible split contract, and this family specializes to that
//! system at comptime.

const std = @import("std");
const testing = std.testing;

pub fn Leapfrog(comptime System: type) type {
    comptime validateSystem(System);

    const SystemDeclType = switch (@typeInfo(System)) {
        .pointer => |pointer| pointer.child,
        else => System,
    };

    return struct {
        pub fn advance(
            allocator: std.mem.Allocator,
            system: System,
            dt: f64,
        ) !void {
            try SystemDeclType.Leapfrog.first(allocator, system, dt);
            try SystemDeclType.Leapfrog.second(allocator, system, dt);
            if (@hasDecl(SystemDeclType.Leapfrog, "applyBoundary")) {
                SystemDeclType.Leapfrog.applyBoundary(system);
            }
        }
    };
}

fn validateSystem(comptime System: type) void {
    const SystemDeclType = switch (@typeInfo(System)) {
        .pointer => |pointer| pointer.child,
        else => System,
    };

    if (!@hasDecl(SystemDeclType, "Leapfrog")) {
        @compileError("Leapfrog requires the system to declare `pub const Leapfrog = struct { ... }`");
    }

    const Split = SystemDeclType.Leapfrog;
    if (@TypeOf(Split) != type) {
        @compileError("Leapfrog: system `Leapfrog` declaration must be a type");
    }

    inline for (.{ "first", "second" }) |name| {
        if (!@hasDecl(Split, name)) {
            @compileError("Leapfrog requires `pub fn " ++ name ++ "(std.mem.Allocator, System, f64) !void`");
        }

        const info = @typeInfo(@TypeOf(@field(Split, name)));
        if (info != .@"fn") {
            @compileError("Leapfrog: `" ++ name ++ "` must be a function");
        }
        const fn_info = info.@"fn";

        if (fn_info.params.len != 3) {
            @compileError("Leapfrog: `" ++ name ++ "` must take exactly 3 parameters (std.mem.Allocator, System, f64)");
        }
        if (fn_info.params[0].type != std.mem.Allocator) {
            @compileError("Leapfrog: `" ++ name ++ "` parameter 0 must be std.mem.Allocator");
        }
        if (fn_info.params[1].type != System) {
            @compileError("Leapfrog: `" ++ name ++ "` parameter 1 must match the chosen system handle type");
        }
        if (fn_info.params[2].type != f64) {
            @compileError("Leapfrog: `" ++ name ++ "` parameter 2 must be f64");
        }

        const ret = fn_info.return_type orelse
            @compileError("Leapfrog: `" ++ name ++ "` must have a known return type");
        const ret_info = @typeInfo(ret);
        if (ret_info != .error_union or ret_info.error_union.payload != void) {
            @compileError("Leapfrog: `" ++ name ++ "` must return !void");
        }
    }

    if (@hasDecl(Split, "applyBoundary")) {
        const info = @typeInfo(@TypeOf(Split.applyBoundary));
        if (info != .@"fn") {
            @compileError("Leapfrog: `applyBoundary` must be a function");
        }
        const fn_info = info.@"fn";
        if (fn_info.params.len != 1 or fn_info.params[0].type != System) {
            @compileError("Leapfrog: `applyBoundary` must take exactly one system-handle parameter");
        }
        if ((fn_info.return_type orelse void) != void) {
            @compileError("Leapfrog: `applyBoundary` must return void");
        }
    }
}

const MockSystem = struct {
    position: f64,
    velocity: f64,

    pub const Leapfrog = struct {
        pub fn first(_: std.mem.Allocator, system: *MockSystem, dt: f64) !void {
            system.position += dt * system.velocity;
        }

        pub fn second(_: std.mem.Allocator, system: *MockSystem, dt: f64) !void {
            system.velocity += dt * (-system.position);
        }
    };
};

test "Leapfrog accepts a conforming system contract" {
    _ = Leapfrog(*MockSystem);
}

test "Leapfrog composes first then second split steps" {
    const Method = Leapfrog(*MockSystem);
    var system = MockSystem{ .position = 0.0, .velocity = 1.0 };

    try Method.advance(testing.allocator, &system, 0.1);
    try testing.expectApproxEqAbs(0.1, system.position, 1e-15);
    try testing.expectApproxEqAbs(0.99, system.velocity, 1e-15);
}

test "Leapfrog preserves harmonic-oscillator energy up to bounded oscillation" {
    const Method = Leapfrog(*MockSystem);
    var system = MockSystem{ .position = 1.0, .velocity = 0.0 };

    const initial_energy = 0.5 * (system.position * system.position + system.velocity * system.velocity);
    const dt = 0.01;
    for (0..10_000) |_| {
        try Method.advance(testing.allocator, &system, dt);
    }

    const final_energy = 0.5 * (system.position * system.position + system.velocity * system.velocity);
    try testing.expect(@abs(final_energy - initial_energy) < 1e-2);
}
