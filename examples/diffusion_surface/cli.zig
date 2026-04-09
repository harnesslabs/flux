//! `flux-examples diffusion-surface` subcommand: heat equation on a sphere.

const std = @import("std");
const common = @import("examples_common");
const surface = @import("example_physics");

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = surface.Config{};
    var co = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

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

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    _ = try surface.run(allocator, config, stderr);
}

fn applyCommon(cfg: *surface.Config, co: common.Common) void {
    common.applySharedFields(cfg, co);
    // diffusion_surface derives dt from final_time / steps; --dt is treated as
    // an override that pins final_time = dt * steps.
    if (co.dt) |dt_value| {
        cfg.final_time = dt_value * @as(f64, @floatFromInt(cfg.steps));
    }
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples diffusion-surface — heat equation on a sphere
        \\
        \\  usage: flux-examples diffusion-surface [options]
        \\
        \\  options:
        \\    --refinement N    icosahedral refinement level (default: 0)
        \\    --steps N         backward-Euler steps (default: 8)
        \\    --dt-scale C      parabolic scale hint (default: 0.1)
        \\    --dt DT           explicit timestep — pins final-time = dt * steps
        \\    --final-time T    total simulated time (default: 0.05)
        \\    --output DIR      output directory for VTK files
        \\    --frames N        number of snapshots to write (default: 4)
        \\
    , .{});
}
