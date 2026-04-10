//! `flux-examples euler-3d` subcommand: 3D incompressible Euler.

const std = @import("std");
const common = @import("examples_common");
const euler = @import("example_physics");

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = euler.Config{};
    var co = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (try common.tryBox3Flag(&parser, arg, &config)) continue;
        if (try parser.tryCommon(arg, &co)) continue;
        std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
        printUsage();
        std.process.exit(2);
    }

    if (co.help) {
        printUsage();
        return;
    }

    applyCommon(&config, co);

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try euler.run(allocator, config, stderr);
    try stderr.print(
        "elapsed={d:.3}s snapshots={d} drift={e}\n",
        .{
            result.elapsed_s,
            result.snapshot_count,
            @abs(result.helicity_final - result.helicity_initial),
        },
    );
}

fn applyCommon(cfg: *euler.Config, co: common.Common) void {
    common.applySharedFields(cfg, co);
    if (co.dt) |v| cfg.dt = v;
    if (co.frames) |frames| {
        cfg.output_interval = common.framesToInterval(cfg.steps, frames);
    }
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples euler-3d — 3D incompressible Euler with helicity
        \\
        \\  usage: flux-examples euler-3d [options]
        \\
        \\  options:
        \\    --steps N             number of timesteps (default: 1000)
        \\    --nx N                x-axis cube count (default: 2)
        \\    --ny N                y-axis cube count (default: 2)
        \\    --nz N                z-axis cube count (default: 2)
        \\    --width L             domain width (default: 1.0)
        \\    --height L            domain height (default: 1.0)
        \\    --depth L             domain depth (default: 1.0)
        \\    --dt DT               timestep size (default: 0.01)
        \\    --output DIR          output directory for VTK snapshots
        \\    --output-interval N   snapshot cadence in steps (overrides --frames)
        \\    --frames N            number of evenly-spaced snapshots
        \\
    , .{});
}
