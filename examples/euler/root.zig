const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const system_mod = @import("system.zig");

pub const system_api = system_mod;

const SnapshotCadence = union(enum) {
    disabled,
    interval: u32,
    frames: u32,
};

pub fn Config(comptime dim: u8) type {
    return struct {
        steps: u32 = 1000,
        counts: [dim]u32 = defaultCounts(dim),
        extents: [dim]f64 = @splat(1.0),
        courant: f64 = 0.1,
        time_step_override: ?f64 = null,
        output_dir: ?[]const u8 = null,
        snapshot_cadence: SnapshotCadence = .{ .frames = 50 },

        pub fn timeStep(self: @This()) f64 {
            if (self.time_step_override) |value| return value;
            return self.courant * minSpacing(dim, self.counts, self.extents);
        }

        pub fn snapshotInterval(self: @This()) ?u32 {
            if (self.output_dir == null) return null;
            return switch (self.snapshot_cadence) {
                .disabled => null,
                .interval => |value| value,
                .frames => |value| @max(@as(u32, 1), common.framesToInterval(self.steps, value)),
            };
        }
    };
}

pub fn RunResult(comptime Summary: type) type {
    return struct {
        elapsed_s: f64,
        snapshot_count: u32,
        summary: Summary,
    };
}

pub const ConservationSummary = struct {
    invariant_initial: f64,
    invariant_final: f64,
};

pub fn runGaussian(
    allocator: std.mem.Allocator,
    config: Config(2),
    writer: *std.Io.Writer,
) !RunResult(ConservationSummary) {
    return try runScenario(2, .gaussian, allocator, config, writer);
}

pub fn runDipole(
    allocator: std.mem.Allocator,
    config: Config(2),
    writer: *std.Io.Writer,
) !RunResult(ConservationSummary) {
    return try runScenario(2, .dipole, allocator, config, writer);
}

pub fn runReference(
    allocator: std.mem.Allocator,
    config: Config(3),
    writer: *std.Io.Writer,
) !RunResult(ConservationSummary) {
    return try runScenario(3, .reference, allocator, config, writer);
}

const Scenario = enum {
    gaussian,
    dipole,
    reference,
};

fn runScenario(
    comptime dim: u8,
    comptime scenario: Scenario,
    allocator: std.mem.Allocator,
    config: Config(dim),
    writer: *std.Io.Writer,
) !RunResult(ConservationSummary) {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Euler = system_mod.Euler(dim, Mesh);

    var mesh = try Mesh.cartesian(allocator, config.counts, config.extents);
    defer mesh.deinit(allocator);

    var system = switch (scenario) {
        .gaussian => try Euler.gaussian(allocator, &mesh),
        .dipole => try Euler.dipole(allocator, &mesh),
        .reference => try Euler.reference(allocator, &mesh),
    };
    defer system.deinit(allocator);

    const initial = try measureInvariant(dim, allocator, &system);
    const result = if (config.snapshotInterval()) |interval| blk: {
        if (dim == 2) {
            var evolution = try flux.evolution.Evolution(*Euler, Euler.Explicit).config()
                .dt(config.timeStep())
                .steps(config.steps)
                .listen(flux.listeners.Progress(writer))
                .listen(
                    flux.listeners.Snapshots(*Euler)
                        .field(.stream_function)
                        .field(.vorticity)
                        .field(.tracer)
                        .field(.velocity)
                        .measurement(.circulation)
                        .directory(config.output_dir.?)
                        .baseName(baseName(scenario))
                        .everySteps(interval),
                )
                .init(allocator, &system);
            defer evolution.deinit();

            const run_result = try evolution.run();
            break :blk RunResult(ConservationSummary){
                .elapsed_s = run_result.elapsed_s,
                .snapshot_count = snapshotCount(config.steps, interval),
                .summary = .{
                    .invariant_initial = initial,
                    .invariant_final = try measureInvariant(dim, allocator, &system),
                },
            };
        }

        var evolution = try flux.evolution.Evolution(*Euler, Euler.Explicit).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .listen(
                flux.listeners.Snapshots(*Euler)
                    .field(.velocity)
                    .field(.vorticity)
                    .measurement(.helicity)
                    .directory(config.output_dir.?)
                    .baseName(baseName(scenario))
                    .everySteps(interval),
            )
            .init(allocator, &system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(ConservationSummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = snapshotCount(config.steps, interval),
            .summary = .{
                .invariant_initial = initial,
                .invariant_final = try measureInvariant(dim, allocator, &system),
            },
        };
    } else blk: {
        var evolution = try flux.evolution.Evolution(*Euler, Euler.Explicit).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .init(allocator, &system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(ConservationSummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = 0,
            .summary = .{
                .invariant_initial = initial,
                .invariant_final = try measureInvariant(dim, allocator, &system),
            },
        };
    };
    return result;
}

fn measureInvariant(comptime dim: u8, allocator: std.mem.Allocator, system: anytype) !f64 {
    if (dim == 2) return system.measurement(allocator, .circulation);
    return system.measurement(allocator, .helicity);
}

fn baseName(comptime scenario: Scenario) []const u8 {
    return switch (scenario) {
        .gaussian => "euler_gaussian",
        .dipole => "euler_dipole",
        .reference => "euler_reference_3d",
    };
}

fn defaultCounts(comptime dim: u8) [dim]u32 {
    return switch (dim) {
        2 => .{ 16, 16 },
        3 => .{ 2, 2, 2 },
        else => @compileError("euler examples only support dimensions 2 and 3"),
    };
}

fn minSpacing(comptime dim: u8, counts: [dim]u32, extents: [dim]f64) f64 {
    var h = extents[0] / @as(f64, @floatFromInt(counts[0]));
    inline for (1..dim) |axis| {
        h = @min(h, extents[axis] / @as(f64, @floatFromInt(counts[axis])));
    }
    return h;
}

fn snapshotCount(steps: u32, interval: u32) u32 {
    std.debug.assert(interval > 0);
    const trailing: u32 = if (steps % interval == 0) 0 else 1;
    return 1 + (steps / interval) + trailing;
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
