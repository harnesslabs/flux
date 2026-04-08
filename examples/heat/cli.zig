//! `flux-examples heat` subcommand: implicit heat equation on a unit square.

const std = @import("std");
const common = @import("examples_common");
const heat = @import("heat.zig");

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = heat.Config{};
    var co = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next(); // discard subcommand name

    while (parser.next()) |arg| {
        if (try parser.tryCommon(arg, &co)) continue;
        if (eql(arg, "--dt-scale")) {
            config.dt_scale = try parser.parseF64("--dt-scale");
        } else {
            std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
            printUsage();
            std.process.exit(2);
        }
    }

    if (co.help) {
        printUsage();
        return;
    }

    applyCommon(&config, co);

    const stderr = (std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter();
    const result = try heat.run(allocator, config, stderr);
    try stderr.print(
        "elapsed={d:.3}s steps={d} snapshots={d} l2_error={e} output={s}\n",
        .{ result.elapsed_s, result.steps, result.snapshot_count, result.l2_error, config.output_dir },
    );
}

fn applyCommon(cfg: *heat.Config, co: common.Common) void {
    if (co.steps) |v| cfg.steps = v;
    if (co.grid) |v| cfg.grid = v;
    if (co.domain) |v| cfg.domain = v;
    if (co.frames) |v| cfg.frames = v;
    if (co.output_dir) |v| cfg.output_dir = v;
    if (co.dt) |v| cfg.dt_override = v;
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples heat — implicit heat equation on the unit square
        \\
        \\  usage: flux-examples heat [options]
        \\
        \\  options:
        \\    --grid N          grid cells per side (default: 8)
        \\    --steps N         backward-Euler steps (default: 8)
        \\    --domain L        square domain side length (default: 1.0)
        \\    --dt-scale C      parabolic timestep scale dt = C·h² (default: 0.1)
        \\    --dt DT           explicit timestep override (replaces dt-scale derived dt)
        \\    --output DIR      output directory for VTK files (default: output/heat)
        \\    --frames N        number of snapshots to write (default: 4)
        \\
    , .{});
}
