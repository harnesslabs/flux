const std = @import("std");
const common = @import("examples_common");

pub const RunResult = struct {
    elapsed_s: f64,
    steps: u32,
    snapshot_count: u32,
    l2_error: f64,
};

pub const RunConfig = struct {
    steps: u32,
    dt: f64,
    final_time: f64,
    frames: u32,
    output_dir: []const u8,
    output_base_name: []const u8,
    progress_writer: ?*std.Io.Writer = null,
};

pub fn runTimeSteppedDiffusion(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    state_values: []f64,
    config: RunConfig,
    exact_initializer: anytype,
    stepper_builder: anytype,
    error_measure: anytype,
) !RunResult {
    std.debug.assert(config.steps > 0);
    std.debug.assert(config.dt > 0.0);
    std.debug.assert(config.final_time > 0.0);

    const scratch_len = stepper_builder.scratchLen();
    const exact = try allocator.alloc(f64, state_values.len);
    defer allocator.free(exact);
    const rhs = try allocator.alloc(f64, scratch_len);
    defer allocator.free(rhs);
    const solution = try allocator.alloc(f64, scratch_len);
    defer allocator.free(solution);

    stepper_builder.seedSolution(state_values, solution);
    const stepper = stepper_builder.makeStepper(state_values, rhs, solution);

    const loop_result = try common.runExactFieldLoop(
        MeshType,
        allocator,
        mesh,
        state_values,
        exact,
        .{
            .steps = config.steps,
            .dt = config.dt,
            .final_time = config.final_time,
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = config.output_base_name,
            .progress_writer = config.progress_writer,
        },
        exact_initializer,
        stepper,
    );

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .steps = config.steps,
        .snapshot_count = loop_result.snapshot_count,
        .l2_error = error_measure.compute(mesh, state_values, exact),
    };
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
