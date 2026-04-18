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

pub fn CavityConfig(comptime dim: u8) type {
    return struct {
        steps: u32 = 1000,
        counts: [dim]u32 = defaultCounts(dim),
        extents: [dim]f64 = @splat(1.0),
        courant: f64 = 0.1,
        time_step_override: ?f64 = null,
        output_dir: ?[]const u8 = null,
        snapshot_cadence: SnapshotCadence = defaultSnapshotCadence(dim),
        reference: bool = false,
        boundary: system_mod.BoundaryCondition = .pec,

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

        pub fn cavityOptions(self: @This()) system_mod.CavityOptions(dim) {
            return .{
                .extents = self.extents,
                .time_step = self.timeStep(),
                .boundary = self.boundary,
            };
        }

        pub fn measurementProvider(self: @This()) system_mod.CavityMeasurementProvider(dim) {
            return .{
                .extents = self.extents,
                .time_step = self.timeStep(),
            };
        }

        pub fn snapshotBaseName(self: @This()) []const u8 {
            return if (self.reference) "cavity_reference" else "cavity";
        }
    };
}

pub const DipoleConfig2D = struct {
    steps: u32 = 1000,
    counts: [2]u32 = .{ 32, 32 },
    extents: [2]f64 = .{ 1.0, 1.0 },
    courant: f64 = 0.1,
    time_step_override: ?f64 = null,
    frequency_hz: f64 = 0.0,
    amplitude: f64 = 1.0,
    output_dir: ?[]const u8 = null,
    snapshot_cadence: SnapshotCadence = .{ .frames = 100 },
    boundary: system_mod.BoundaryCondition = .pec,

    pub fn timeStep(self: @This()) f64 {
        if (self.time_step_override) |value| return value;
        return self.courant * minSpacing(2, self.counts, self.extents);
    }

    pub fn sourceFrequency(self: @This()) f64 {
        if (self.frequency_hz != 0.0) return self.frequency_hz;
        return 1.0 / (2.0 * self.extents[0]);
    }

    pub fn snapshotInterval(self: @This()) ?u32 {
        if (self.output_dir == null) return null;
        return switch (self.snapshot_cadence) {
            .disabled => null,
            .interval => |value| value,
            .frames => |value| @max(@as(u32, 1), common.framesToInterval(self.steps, value)),
        };
    }

    pub fn dipoleOptions(self: @This()) system_mod.DipoleOptions(2) {
        return .{
            .center = .{ 0.5 * self.extents[0], 0.5 * self.extents[1] },
            .frequency_hz = self.sourceFrequency(),
            .amplitude = self.amplitude,
            .boundary = self.boundary,
        };
    }
};

pub fn RunResult(comptime Summary: type) type {
    return struct {
        elapsed_s: f64,
        snapshot_count: u32,
        summary: Summary,
    };
}

pub const CavitySummary = struct {
    energy_final: f64,
    electric_l2_final: ?f64 = null,
    magnetic_l2_final: ?f64 = null,
};

pub fn runDipole(
    allocator: std.mem.Allocator,
    config: DipoleConfig2D,
    writer: *std.Io.Writer,
) !RunResult(CavitySummary) {
    const Mesh = flux.topology.Mesh(2, 2);
    const Maxwell = system_mod.Maxwell(2, Mesh);

    var mesh = try Mesh.cartesian(allocator, config.counts, config.extents);
    defer mesh.deinit(allocator);

    var system = try Maxwell.dipole(allocator, &mesh, config.dipoleOptions());
    defer system.deinit(allocator);

    const result = if (config.snapshotInterval()) |interval| blk: {
        var evolution = try flux.evolution.Evolution(*Maxwell, flux.integrators.Leapfrog).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .listen(
                flux.listeners.Snapshots(*Maxwell)
                    .field(.electric)
                    .field(.magnetic)
                    .measurement(.energy)
                    .directory(config.output_dir.?)
                    .baseName("dipole")
                    .everySteps(interval),
            )
            .init(allocator, &system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(CavitySummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = snapshotCount(config.steps, interval),
            .summary = .{
                .energy_final = try system.measurement(allocator, .energy),
            },
        };
    } else blk: {
        var evolution = try flux.evolution.Evolution(*Maxwell, flux.integrators.Leapfrog).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .init(allocator, &system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(CavitySummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = 0,
            .summary = .{
                .energy_final = try system.measurement(allocator, .energy),
            },
        };
    };
    return result;
}

pub fn runCavity(
    comptime dim: u8,
    allocator: std.mem.Allocator,
    config: CavityConfig(dim),
    writer: *std.Io.Writer,
) !RunResult(CavitySummary) {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Maxwell = system_mod.Maxwell(dim, Mesh);

    var mesh = try Mesh.cartesian(allocator, config.counts, config.extents);
    defer mesh.deinit(allocator);

    var system = try Maxwell.cavity(allocator, &mesh, config.cavityOptions());
    defer system.deinit(allocator);

    return if (config.reference)
        try runCavityMode(true, allocator, config, writer, Maxwell, &system)
    else
        try runCavityMode(false, allocator, config, writer, Maxwell, &system);
}

fn runCavityMode(
    comptime with_reference: bool,
    allocator: std.mem.Allocator,
    config: anytype,
    writer: *std.Io.Writer,
    comptime Maxwell: type,
    system: *Maxwell,
) !RunResult(CavitySummary) {
    const provider = if (with_reference) config.measurementProvider() else {};
    const result = if (config.snapshotInterval()) |interval| blk: {
        var evolution_config = flux.evolution.Evolution(*Maxwell, flux.integrators.Leapfrog).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer));
        if (with_reference) {
            evolution_config = evolution_config.listen(
                flux.listeners.Snapshots(*Maxwell)
                    .measureWith(provider)
                    .field(.electric)
                    .field(.magnetic)
                    .measurement(.energy)
                    .measurement(.electric_l2)
                    .measurement(.magnetic_l2)
                    .directory(config.output_dir.?)
                    .baseName(config.snapshotBaseName())
                    .everySteps(interval),
            );
        } else {
            evolution_config = evolution_config.listen(
                flux.listeners.Snapshots(*Maxwell)
                    .field(.electric)
                    .field(.magnetic)
                    .measurement(.energy)
                    .directory(config.output_dir.?)
                    .baseName(config.snapshotBaseName())
                    .everySteps(interval),
            );
        }

        var evolution = try evolution_config.init(allocator, system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(CavitySummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = snapshotCount(config.steps, interval),
            .summary = try cavitySummary(with_reference, allocator, system, evolution.currentTime(), provider),
        };
    } else blk: {
        var evolution = try flux.evolution.Evolution(*Maxwell, flux.integrators.Leapfrog).config()
            .dt(config.timeStep())
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .init(allocator, system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(CavitySummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = 0,
            .summary = try cavitySummary(with_reference, allocator, system, evolution.currentTime(), provider),
        };
    };
    return result;
}

fn cavitySummary(
    comptime with_reference: bool,
    allocator: std.mem.Allocator,
    system: anytype,
    time: f64,
    provider: anytype,
) !CavitySummary {
    if (with_reference) {
        return .{
            .energy_final = try provider.measurement(allocator, system, .energy, time),
            .electric_l2_final = try provider.measurement(allocator, system, .electric_l2, time),
            .magnetic_l2_final = try provider.measurement(allocator, system, .magnetic_l2, time),
        };
    }

    return .{
        .energy_final = try system.measurement(allocator, .energy),
    };
}

fn defaultCounts(comptime dim: u8) [dim]u32 {
    return switch (dim) {
        2 => .{ 32, 32 },
        3 => .{ 2, 2, 2 },
        else => @compileError("new_maxwell only supports dimensions 2 and 3"),
    };
}

fn defaultSnapshotCadence(comptime dim: u8) SnapshotCadence {
    return switch (dim) {
        2 => .{ .frames = 100 },
        3 => .disabled,
        else => @compileError("new_maxwell only supports dimensions 2 and 3"),
    };
}

fn minSpacing(comptime dim: u8, counts: [dim]u32, extents: [dim]f64) f64 {
    var spacing = extents[0] / @as(f64, @floatFromInt(counts[0]));
    inline for (1..dim) |axis| {
        spacing = @min(spacing, extents[axis] / @as(f64, @floatFromInt(counts[axis])));
    }
    return spacing;
}

fn snapshotCount(steps: u32, interval: u32) u32 {
    std.debug.assert(interval > 0);
    const trailing: u32 = if (steps % interval == 0) 0 else 1;
    return 1 + (steps / interval) + trailing;
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
