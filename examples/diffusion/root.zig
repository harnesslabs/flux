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

pub const PlaneConfig = struct {
    steps: u32 = 8,
    counts: [2]u32 = .{ 8, 8 },
    extents: [2]f64 = .{ 1.0, 1.0 },
    dt_scale: f64 = 0.1,
    time_step_override: ?f64 = null,
    output_dir: ?[]const u8 = null,
    snapshot_cadence: SnapshotCadence = .{ .frames = 4 },

    pub fn timeStep(self: @This()) f64 {
        if (self.time_step_override) |value| return value;
        const h = @min(
            self.extents[0] / @as(f64, @floatFromInt(self.counts[0])),
            self.extents[1] / @as(f64, @floatFromInt(self.counts[1])),
        );
        return self.dt_scale * h * h;
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

pub const SphereConfig = struct {
    steps: u32 = 8,
    refinement: u32 = 0,
    radius: f64 = 1.0,
    final_time: f64 = 0.05,
    output_dir: ?[]const u8 = null,
    snapshot_cadence: SnapshotCadence = .{ .frames = 4 },

    pub fn timeStep(self: @This()) f64 {
        return self.final_time / @as(f64, @floatFromInt(self.steps));
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

pub fn RunResult(comptime Summary: type) type {
    return struct {
        elapsed_s: f64,
        snapshot_count: u32,
        summary: Summary,
    };
}

pub const DiffusionSummary = struct {
    l2_error_final: f64,
};

pub fn runPlane(
    allocator: std.mem.Allocator,
    config: PlaneConfig,
    writer: *std.Io.Writer,
) !RunResult(DiffusionSummary) {
    const Mesh = flux.topology.Mesh(2, 2);
    const Diffusion = system_mod.Diffusion(.plane, Mesh);

    var mesh = try Mesh.cartesian(allocator, config.counts, config.extents);
    defer mesh.deinit(allocator);

    var system = try Diffusion.plane(allocator, &mesh, config.timeStep());
    defer system.deinit(allocator);

    const provider = system_mod.ReferenceMeasurementProvider(.plane, Mesh){};
    return try runMode(
        allocator,
        writer,
        config.steps,
        config.timeStep(),
        config.snapshotInterval(),
        config.output_dir,
        "diffusion_plane",
        Diffusion,
        &system,
        provider,
    );
}

pub fn runSphere(
    allocator: std.mem.Allocator,
    config: SphereConfig,
    writer: *std.Io.Writer,
) !RunResult(DiffusionSummary) {
    const Mesh = flux.topology.Mesh(3, 2);
    const Diffusion = system_mod.Diffusion(.sphere, Mesh);

    var mesh = try Mesh.sphere(allocator, config.radius, config.refinement);
    defer mesh.deinit(allocator);

    var system = try Diffusion.sphere(allocator, &mesh, config.timeStep());
    defer system.deinit(allocator);

    const provider = system_mod.ReferenceMeasurementProvider(.sphere, Mesh){};
    return try runMode(
        allocator,
        writer,
        config.steps,
        config.timeStep(),
        config.snapshotInterval(),
        config.output_dir,
        "diffusion_sphere",
        Diffusion,
        &system,
        provider,
    );
}

fn runMode(
    allocator: std.mem.Allocator,
    writer: *std.Io.Writer,
    steps: u32,
    time_step: f64,
    snapshot_interval: ?u32,
    output_dir: ?[]const u8,
    base_name: []const u8,
    comptime Diffusion: type,
    system: *Diffusion,
    provider: anytype,
) !RunResult(DiffusionSummary) {
    const result = if (snapshot_interval) |interval| blk: {
        var evolution = try flux.evolution.Evolution(*Diffusion, Diffusion.BackwardEuler).config()
            .dt(time_step)
            .steps(steps)
            .listen(flux.listeners.Progress(writer))
            .listen(
                flux.listeners.Snapshots(*Diffusion)
                    .measureWith(provider)
                    .field(.temperature)
                    .measurement(.l2_error)
                    .directory(output_dir.?)
                    .baseName(base_name)
                    .everySteps(interval),
            )
            .init(allocator, system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(DiffusionSummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = snapshotCount(steps, interval),
            .summary = .{
                .l2_error_final = try provider.measurement(allocator, system, .l2_error, evolution.currentTime()),
            },
        };
    } else blk: {
        var evolution = try flux.evolution.Evolution(*Diffusion, Diffusion.BackwardEuler).config()
            .dt(time_step)
            .steps(steps)
            .listen(flux.listeners.Progress(writer))
            .init(allocator, system);
        defer evolution.deinit();

        const run_result = try evolution.run();
        break :blk RunResult(DiffusionSummary){
            .elapsed_s = run_result.elapsed_s,
            .snapshot_count = 0,
            .summary = .{
                .l2_error_final = try provider.measurement(allocator, system, .l2_error, evolution.currentTime()),
            },
        };
    };
    return result;
}

fn snapshotCount(steps: u32, interval: u32) u32 {
    std.debug.assert(interval > 0);
    const trailing: u32 = if (steps % interval == 0) 0 else 1;
    return 1 + (steps / interval) + trailing;
}

fn planeConvergenceConfig(count: u32) PlaneConfig {
    return .{
        .steps = count * count,
        .counts = .{ count, count },
        .time_step_override = 0.1 / (@as(f64, @floatFromInt(count)) * @as(f64, @floatFromInt(count))),
    };
}

test "plane diffusion reaches second-order spatial convergence under refinement" {
    const allocator = std.testing.allocator;
    var coarse_writer = std.Io.Writer.Allocating.init(allocator);
    defer coarse_writer.deinit();
    var fine_writer = std.Io.Writer.Allocating.init(allocator);
    defer fine_writer.deinit();
    var finer_writer = std.Io.Writer.Allocating.init(allocator);
    defer finer_writer.deinit();

    const coarse = try runPlane(allocator, planeConvergenceConfig(8), &coarse_writer.writer);
    const fine = try runPlane(allocator, planeConvergenceConfig(16), &fine_writer.writer);
    const finer = try runPlane(allocator, planeConvergenceConfig(32), &finer_writer.writer);

    const errors = [_]f64{
        coarse.summary.l2_error_final,
        fine.summary.l2_error_final,
        finer.summary.l2_error_final,
    };
    const rates = try flux.evolution.reference.empiricalRates(allocator, &errors, 2.0);
    defer allocator.free(rates);

    try std.testing.expect(fine.summary.l2_error_final < coarse.summary.l2_error_final);
    try std.testing.expect(finer.summary.l2_error_final < fine.summary.l2_error_final);
    for (rates) |rate| {
        try std.testing.expect(rate > 1.75);
    }
}

test "sphere diffusion error decreases under refinement" {
    const allocator = std.testing.allocator;
    var coarse_writer = std.Io.Writer.Allocating.init(allocator);
    defer coarse_writer.deinit();
    var fine_writer = std.Io.Writer.Allocating.init(allocator);
    defer fine_writer.deinit();

    const coarse = try runSphere(allocator, .{ .refinement = 1 }, &coarse_writer.writer);
    const fine = try runSphere(allocator, .{ .refinement = 2 }, &fine_writer.writer);

    try std.testing.expect(fine.summary.l2_error_final < coarse.summary.l2_error_final);
    try std.testing.expect(fine.summary.l2_error_final < 5e-3);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
