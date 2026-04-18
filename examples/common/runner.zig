//! Shared run-loop helpers for examples with optional reference-study capture.
//!
//! This owns the time loop, snapshot cadence, reference-field visualization, and
//! convergence-study bookkeeping. Mesh assembly, state evolution, and physics
//! invariants stay in the family modules.

const std = @import("std");
const flux = @import("flux");
const snapshot = @import("snapshot.zig");
const cli = @import("cli.zig");
const progress_mod = @import("progress.zig");

const flux_io = flux.io;

pub const RunLoopConfig = struct {
    steps: u32,
    final_time: f64,
    frames: u32,
    output_dir: []const u8,
    output_base_name: []const u8,
    capture_initial: bool = true,
    emit_final: bool = true,
    plan_override: ?snapshot.Plan = null,
    progress_writer: ?*std.Io.Writer = null,
};

pub const RunLoopResult = struct {
    elapsed_s: f64,
    snapshot_count: u32,
};

pub fn ReferenceFieldRenderer(comptime MeshType: type) type {
    return struct {
        mesh: *const MeshType,
        state: []const f64,
        exact: []const f64,

        pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
            const error_values = try allocator.alloc(f64, self.state.len);
            defer allocator.free(error_values);

            for (self.state, self.exact, error_values) |approx_value, exact_value, *error_value| {
                error_value.* = approx_value - exact_value;
            }

            const point_data = [_]flux_io.DataArraySlice{
                .{ .name = "temperature", .values = self.state },
                .{ .name = "temperature_exact", .values = self.exact },
                .{ .name = "temperature_error", .values = error_values },
            };
            try flux_io.write(writer, self.mesh.*, &point_data, &.{});
        }
    };
}

fn buildPlan(config: RunLoopConfig) snapshot.Plan {
    if (config.plan_override) |plan| return plan;
    if (config.frames > 0) {
        const interval = cli.framesToInterval(config.steps, config.frames);
        return snapshot.Plan.fromInterval(config.steps, interval, .{
            .emit_final = config.emit_final,
        });
    }
    return snapshot.Plan.disabled();
}

fn runCapturedLoop(
    allocator: std.mem.Allocator,
    evolution: anytype,
    config: RunLoopConfig,
    capturer: anytype,
) !RunLoopResult {
    const EvolutionNs = @TypeOf(evolution.*).EvolutionNamespace;
    const plan = buildPlan(config);

    std.debug.assert(evolution.configuredSteps() == config.steps);
    const expected_final_time = evolution.currentTime() +
        evolution.timeStep() * @as(f64, @floatFromInt(config.steps));
    std.debug.assert(
        std.math.approxEqAbs(f64, expected_final_time, config.final_time, 1e-12) or
            std.math.approxEqRel(f64, expected_final_time, config.final_time, 1e-12),
    );

    var series = try snapshot.Series.init(
        allocator,
        config.output_dir,
        config.output_base_name,
        plan,
    );
    defer series.deinit();

    const SnapshotListener = struct {
        series_ptr: *snapshot.Series,
        inner: @TypeOf(capturer),
        capture_initial: bool,
        total_steps: u32,

        pub fn onRunBegin(self: *@This(), event: EvolutionNs.Event) !void {
            if (!self.capture_initial or !self.series_ptr.enabled()) return;
            _ = event;
            try self.inner.capture(self.series_ptr, 0.0);
        }

        pub fn onStep(self: *@This(), event: EvolutionNs.Event) !void {
            if (!self.series_ptr.enabled()) return;
            if (!self.series_ptr.dueAt(event.step_index, self.total_steps)) return;
            try self.inner.capture(self.series_ptr, event.time);
        }

        pub fn onRunEnd(self: *@This(), event: EvolutionNs.Event) !void {
            _ = event;
            try self.series_ptr.finalize();
        }
    };

    var snapshot_listener = SnapshotListener{
        .series_ptr = &series,
        .inner = capturer,
        .capture_initial = config.capture_initial,
        .total_steps = config.steps,
    };

    var progress = if (config.progress_writer) |writer|
        progress_mod.Progress.init(writer, config.steps)
    else
        null;

    const ProgressListener = struct {
        progress_ptr: *progress_mod.Progress,

        pub fn onStep(self: *@This(), event: EvolutionNs.Event) !void {
            self.progress_ptr.update(event.step_index);
        }

        pub fn onRunEnd(self: *@This(), event: EvolutionNs.Event) !void {
            _ = event;
            self.progress_ptr.finish();
        }
    };

    const run_result = if (progress) |*bar| blk: {
        var progress_listener = ProgressListener{ .progress_ptr = bar };
        break :blk try evolution.runWith(.{
            &snapshot_listener,
            &progress_listener,
        });
    } else try evolution.runWith(.{&snapshot_listener});

    return .{
        .elapsed_s = run_result.elapsed_s,
        .snapshot_count = series.count,
    };
}

pub fn runEvolutionLoop(
    allocator: std.mem.Allocator,
    evolution: anytype,
    config: RunLoopConfig,
    renderer: anytype,
) !RunLoopResult {
    const Capturer = struct {
        renderer_value: @TypeOf(renderer),

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            try series.capture(time, self.renderer_value);
        }
    };

    return runCapturedLoop(
        allocator,
        evolution,
        config,
        Capturer{ .renderer_value = renderer },
    );
}

pub fn runReferenceEvolutionLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    evolution: anytype,
    study: anytype,
    config: RunLoopConfig,
) !RunLoopResult {
    const Capturer = struct {
        mesh_ptr: *const MeshType,
        study_value: @TypeOf(study),

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            self.study_value.fillExact(time);
            try series.capture(time, ReferenceFieldRenderer(MeshType){
                .mesh = self.mesh_ptr,
                .state = self.study_value.stateValues(),
                .exact = self.study_value.exactValues(),
            });
        }
    };

    const result = try runCapturedLoop(
        allocator,
        evolution,
        config,
        Capturer{
            .mesh_ptr = mesh,
            .study_value = study,
        },
    );
    study.fillExact(config.final_time);
    return result;
}

pub fn runConvergenceStudy(
    comptime ResultType: type,
    comptime ParamType: type,
    allocator: std.mem.Allocator,
    params: []const ParamType,
    case_runner: anytype,
) ![]ResultType {
    const results = try allocator.alloc(ResultType, params.len);
    errdefer allocator.free(results);
    const RunnerType = @TypeOf(case_runner);

    for (params, 0..) |param, idx| {
        results[idx] = try RunnerType.run(allocator, param);
    }

    return results;
}

test "runConvergenceStudy preserves parameter order" {
    const Result = struct { value: u32 };
    const Runner = struct {
        pub fn run(_: std.mem.Allocator, value: u32) !Result {
            return .{ .value = value };
        }
    };

    const allocator = std.testing.allocator;
    const params = [_]u32{ 3, 1, 4 };
    const results = try runConvergenceStudy(Result, u32, allocator, &params, Runner{});
    defer allocator.free(results);

    try std.testing.expectEqual(@as(u32, 3), results[0].value);
    try std.testing.expectEqual(@as(u32, 1), results[1].value);
    try std.testing.expectEqual(@as(u32, 4), results[2].value);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
