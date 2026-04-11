//! Shared run-loop helpers for examples with optional exact/reference fields.
//!
//! This owns the time loop, snapshot cadence, exact-field visualization, and
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
    dt: f64,
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

pub fn ExactFieldRenderer(comptime MeshType: type) type {
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

pub fn runStepLoop(
    allocator: std.mem.Allocator,
    config: RunLoopConfig,
    stepper: anytype,
    capturer: anytype,
) !RunLoopResult {
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.dt > 0.0);
    std.debug.assert(config.final_time > 0.0);

    const plan = buildPlan(config);

    var series = try snapshot.Series.init(
        allocator,
        config.output_dir,
        config.output_base_name,
        plan,
    );
    defer series.deinit();

    if (config.capture_initial and series.enabled()) {
        try capturer.capture(&series, 0.0);
    }

    var progress = if (config.progress_writer) |writer|
        progress_mod.Progress.init(writer, config.steps)
    else
        null;

    const start_ns = std.time.nanoTimestamp();
    for (0..config.steps) |step_idx| {
        const next_time = config.dt * @as(f64, @floatFromInt(step_idx + 1));
        try stepper.step(allocator);

        if (series.enabled() and series.dueAt(@intCast(step_idx + 1), config.steps)) {
            try capturer.capture(&series, next_time);
        }

        if (progress) |*bar| {
            bar.update(@intCast(step_idx + 1));
        }
    }
    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

    if (progress) |*bar| {
        bar.finish();
    }
    try series.finalize();

    return .{
        .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
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

    return runStepLoop(
        allocator,
        config,
        evolution,
        Capturer{ .renderer_value = renderer },
    );
}

pub fn runSimulationLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    state_values: []const f64,
    exact_values: ?[]f64,
    config: RunLoopConfig,
    exact_initializer: anytype,
    stepper: anytype,
) !RunLoopResult {
    const Capturer = struct {
        mesh_ptr: *const MeshType,
        state: []const f64,
        exact: ?[]f64,
        exact_init: @TypeOf(exact_initializer),

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            if (self.exact) |exact_values_inner| {
                self.exact_init.fill(self.mesh_ptr, exact_values_inner, time);
                try series.capture(time, ExactFieldRenderer(MeshType){
                    .mesh = self.mesh_ptr,
                    .state = self.state,
                    .exact = exact_values_inner,
                });
            }
        }
    };

    const result = try runStepLoop(
        allocator,
        config,
        stepper,
        Capturer{
            .mesh_ptr = mesh,
            .state = state_values,
            .exact = exact_values,
            .exact_init = exact_initializer,
        },
    );

    if (exact_values) |exact| {
        exact_initializer.fill(mesh, exact, config.final_time);
    }
    return result;
}

pub fn runExactEvolutionLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    evolution: anytype,
    config: RunLoopConfig,
) !RunLoopResult {
    const EvolutionType = @TypeOf(evolution.*);

    const Capturer = struct {
        mesh_ptr: *const MeshType,
        evolution_ptr: *EvolutionType,

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            self.evolution_ptr.fillExact(time);
            try series.capture(time, ExactFieldRenderer(MeshType){
                .mesh = self.mesh_ptr,
                .state = self.evolution_ptr.stateValues(),
                .exact = self.evolution_ptr.exactValues(),
            });
        }
    };

    const result = try runStepLoop(
        allocator,
        config,
        evolution,
        Capturer{
            .mesh_ptr = mesh,
            .evolution_ptr = evolution,
        },
    );
    evolution.fillExact(config.final_time);
    return result;
}

pub fn runExactFieldLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    state_values: []const f64,
    exact_values: []f64,
    config: RunLoopConfig,
    exact_initializer: anytype,
    stepper: anytype,
) !RunLoopResult {
    return runSimulationLoop(
        MeshType,
        allocator,
        mesh,
        state_values,
        exact_values,
        config,
        exact_initializer,
        stepper,
    );
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
