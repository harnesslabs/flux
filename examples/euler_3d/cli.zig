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
        if (eql(arg, "--nx")) {
            config.nx = try parser.parseU32("--nx");
            continue;
        }
        if (eql(arg, "--ny")) {
            config.ny = try parser.parseU32("--ny");
            continue;
        }
        if (eql(arg, "--nz")) {
            config.nz = try parser.parseU32("--nz");
            continue;
        }
        if (eql(arg, "--width")) {
            config.width = try parser.parseF64("--width");
            continue;
        }
        if (eql(arg, "--height")) {
            config.height = try parser.parseF64("--height");
            continue;
        }
        if (eql(arg, "--depth")) {
            config.depth = try parser.parseF64("--depth");
            continue;
        }
        if (eql(arg, "--output-interval")) {
            config.output_interval = try parser.parseU32("--output-interval");
            continue;
        }
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
    if (co.steps) |v| cfg.steps = v;
    if (co.dt) |v| cfg.dt = v;
    // euler_3d's Config stores output_dir as ?[]const u8 (nullable means
    // "don't write snapshots"). The shared flag is non-optional so we lift
    // it back into the option here.
    if (co.output_dir) |v| cfg.output_dir = v;
    if (co.frames) |frames| {
        // Translate the shared --frames flag into the example's
        // output_interval semantics so 3D Euler honors the same convention
        // as 2D simulations.
        if (frames == 0) {
            cfg.output_interval = 0;
        } else {
            cfg.output_interval = @max(@as(u32, 1), cfg.steps / frames);
        }
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
