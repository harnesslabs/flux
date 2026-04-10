//! `flux-examples maxwell-3d` subcommand: 3D Maxwell TM_110 cavity.

const std = @import("std");
const common = @import("examples_common");
const maxwell = @import("example_physics");

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = maxwell.Config{};
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

    return maxwell.runDriver(allocator, config);
}

fn applyCommon(cfg: *maxwell.Config, co: common.Common) void {
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
        \\  flux-examples maxwell-3d — TM_110 cavity resonance on tetrahedra
        \\
        \\  usage: flux-examples maxwell-3d [options]
        \\
        \\  options:
        \\    --nx N                tetrahedral cells in x (default: 2)
        \\    --ny N                tetrahedral cells in y (default: 2)
        \\    --nz N                tetrahedral cells in z (default: 2)
        \\    --width L             cavity width  (default: 1.0)
        \\    --height L            cavity height (default: 1.0)
        \\    --depth L             cavity depth  (default: 1.0)
        \\    --steps N             leapfrog steps (default: 1000)
        \\    --dt DT               fixed timestep (default: 0.01)
        \\    --output DIR          write VTK snapshots into DIR
        \\    --output-interval N   write every N steps (overrides --frames)
        \\    --frames N            number of evenly-spaced snapshots
        \\
    , .{});
}
