//! `flux-examples euler-2d` subcommand: 2D incompressible Euler.

const std = @import("std");
const common = @import("examples_common");
const euler = @import("euler.zig");

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = euler.Config{};
    var co = common.Common{};
    var output_explicit = false;
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (eql(arg, "--demo")) {
            const name = try parser.requireValue("--demo");
            if (eql(name, "gaussian")) {
                config.demo = .gaussian;
            } else if (eql(name, "dipole")) {
                config.demo = .dipole;
            } else {
                std.debug.print("error: unknown demo '{s}'. available: gaussian, dipole\n", .{name});
                std.process.exit(2);
            }
            continue;
        }
        if (eql(arg, "--cfl")) {
            config.cfl = try parser.parseF64("--cfl");
            continue;
        }
        if (try parser.tryCommon(arg, &co)) {
            if (co.output_dir != null) output_explicit = true;
            continue;
        }
        std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
        printUsage();
        std.process.exit(2);
    }

    if (co.help) {
        printUsage();
        return;
    }

    applyCommon(&config, co);
    // Match the previous default-output convention: dipole demos write to a
    // dedicated subdirectory unless the user provided --output explicitly.
    if (!output_explicit and config.demo == .dipole) {
        config.output_dir = "output/euler_dipole";
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try euler.run(allocator, config, stderr);
    const drift = @abs(result.circulation_final - result.circulation_initial);
    try stderr.print(
        "elapsed={d:.3}s snapshots={d} drift={e} output={s}\n",
        .{ result.elapsed_s, result.snapshot_count, drift, config.output_dir },
    );
}

fn applyCommon(cfg: *euler.Config, co: common.Common) void {
    if (co.steps) |v| cfg.steps = v;
    if (co.grid) |v| cfg.grid = v;
    if (co.domain) |v| cfg.domain = v;
    if (co.frames) |v| cfg.frames = v;
    if (co.output_dir) |v| cfg.output_dir = v;
    if (co.dt) |dt_value| {
        // 2D Euler derives dt from cfl * (domain / grid). Solve for the cfl
        // that produces the requested dt so the existing dt() method stays
        // the single source of truth.
        const h = cfg.domain / @as(f64, @floatFromInt(cfg.grid));
        cfg.cfl = dt_value / h;
    }
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples euler-2d — 2D incompressible Euler vorticity-stream
        \\
        \\  usage: flux-examples euler-2d [options]
        \\
        \\  options:
        \\    --demo NAME       gaussian (invariant reference) or dipole (dynamic showcase)
        \\    --grid N          grid cells per side (default: 16)
        \\    --steps N         number of timesteps (default: 1000)
        \\    --domain L        square domain side length (default: 1.0)
        \\    --cfl C           timestep scale dt = C·h (default: 0.1)
        \\    --dt DT           explicit timestep override (computes cfl = dt/h)
        \\    --output DIR      output directory for VTK files (default: output/euler_2d)
        \\    --frames N        number of snapshots to write (default: 50)
        \\
    , .{});
}
