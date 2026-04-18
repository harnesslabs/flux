const std = @import("std");
const euler = @import("new_euler");
const parse = @import("parse.zig");

pub const name = "euler";
pub const summary = "Fresh Euler examples over examples/new_euler";

const Scenario = enum {
    gaussian,
    dipole,
    reference,
};

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var parser = parse.Parser.init(args);
    _ = parser.next();

    const scenario_arg = parser.next() orelse {
        printUsage();
        std.process.exit(2);
    };

    if (parse.eql(scenario_arg, "--help") or parse.eql(scenario_arg, "-h") or parse.eql(scenario_arg, "help")) {
        printUsage();
        return;
    }

    const scenario: Scenario = if (parse.eql(scenario_arg, "gaussian"))
        .gaussian
    else if (parse.eql(scenario_arg, "dipole"))
        .dipole
    else if (parse.eql(scenario_arg, "reference"))
        .reference
    else {
        std.debug.print("error: unknown Euler scenario '{s}'\n\n", .{scenario_arg});
        printUsage();
        std.process.exit(2);
    };

    switch (scenario) {
        .gaussian => try run2DScenario(.gaussian, allocator, &parser),
        .dipole => try run2DScenario(.dipole, allocator, &parser),
        .reference => try runReference(allocator, &parser),
    }
}

fn run2DScenario(comptime scenario: Scenario, allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = euler.Config(2){};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            print2DUsage();
            return;
        }
        if (parse.eql(arg, "--steps")) {
            config.steps = try parser.parseU32("--steps");
            continue;
        }
        if (parse.eql(arg, "--nx")) {
            config.counts[0] = try parser.parseU32("--nx");
            continue;
        }
        if (parse.eql(arg, "--ny")) {
            config.counts[1] = try parser.parseU32("--ny");
            continue;
        }
        if (parse.eql(arg, "--width")) {
            config.extents[0] = try parser.parseF64("--width");
            continue;
        }
        if (parse.eql(arg, "--height")) {
            config.extents[1] = try parser.parseF64("--height");
            continue;
        }
        if (parse.eql(arg, "--courant")) {
            config.courant = try parser.parseF64("--courant");
            continue;
        }
        if (parse.eql(arg, "--dt")) {
            config.time_step_override = try parser.parseF64("--dt");
            continue;
        }
        if (parse.eql(arg, "--output")) {
            config.output_dir = try parser.requireValue("--output");
            continue;
        }
        if (parse.eql(arg, "--frames")) {
            config.snapshot_cadence = .{ .frames = try parser.parseU32("--frames") };
            continue;
        }
        if (parse.eql(arg, "--output-interval")) {
            config.snapshot_cadence = .{ .interval = try parser.parseU32("--output-interval") };
            continue;
        }
        failUnknownFlag(arg, print2DUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = switch (scenario) {
        .gaussian => try euler.runGaussian(allocator, config, stderr),
        .dipole => try euler.runDipole(allocator, config, stderr),
        else => unreachable,
    };
    try stderr.print(
        "euler_{s}: elapsed={d:.3}s snapshots={d} invariant={e} -> {e}\n",
        .{
            if (scenario == .gaussian) "gaussian" else "dipole",
            result.elapsed_s,
            result.snapshot_count,
            result.summary.invariant_initial,
            result.summary.invariant_final,
        },
    );
}

fn runReference(allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = euler.Config(3){};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            print3DUsage();
            return;
        }
        if (parse.eql(arg, "--steps")) {
            config.steps = try parser.parseU32("--steps");
            continue;
        }
        if (try consumeCounts(3, parser, arg, &config.counts)) continue;
        if (try consumeExtents(3, parser, arg, &config.extents)) continue;
        if (parse.eql(arg, "--courant")) {
            config.courant = try parser.parseF64("--courant");
            continue;
        }
        if (parse.eql(arg, "--dt")) {
            config.time_step_override = try parser.parseF64("--dt");
            continue;
        }
        if (parse.eql(arg, "--output")) {
            config.output_dir = try parser.requireValue("--output");
            continue;
        }
        if (parse.eql(arg, "--frames")) {
            config.snapshot_cadence = .{ .frames = try parser.parseU32("--frames") };
            continue;
        }
        if (parse.eql(arg, "--output-interval")) {
            config.snapshot_cadence = .{ .interval = try parser.parseU32("--output-interval") };
            continue;
        }
        failUnknownFlag(arg, print3DUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = try euler.runReference(allocator, config, stderr);
    try stderr.print(
        "euler_reference_3d: elapsed={d:.3}s snapshots={d} invariant={e} -> {e}\n",
        .{
            result.elapsed_s,
            result.snapshot_count,
            result.summary.invariant_initial,
            result.summary.invariant_final,
        },
    );
}

fn consumeCounts(comptime dim: u8, parser: *parse.Parser, arg: []const u8, counts: *[dim]u32) !bool {
    inline for (0..dim) |axis| {
        const flag = switch (axis) {
            0 => "--n0",
            1 => "--n1",
            2 => "--n2",
            else => unreachable,
        };
        if (parse.eql(arg, flag)) {
            counts[axis] = try parser.parseU32(flag);
            return true;
        }
    }
    return false;
}

fn consumeExtents(comptime dim: u8, parser: *parse.Parser, arg: []const u8, extents: *[dim]f64) !bool {
    inline for (0..dim) |axis| {
        const flag = switch (axis) {
            0 => "--extent0",
            1 => "--extent1",
            2 => "--extent2",
            else => unreachable,
        };
        if (parse.eql(arg, flag)) {
            extents[axis] = try parser.parseF64(flag);
            return true;
        }
    }
    return false;
}

fn failUnknownFlag(arg: []const u8, usage: *const fn () void) noreturn {
    std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
    usage();
    std.process.exit(2);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli euler <scenario> [options]
        \\
        \\  scenarios:
        \\    gaussian   2D Gaussian vortex
        \\    dipole     2D vortex dipole
        \\    reference  3D reference mode
        \\
        \\  ask a scenario for details:
        \\    zig build run-new-cli -- euler gaussian --help
        \\    zig build run-new-cli -- euler dipole --help
        \\    zig build run-new-cli -- euler reference --help
        \\
    , .{});
}

fn print2DUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli euler <gaussian|dipole> [options]
        \\
        \\  options:
        \\    --steps N
        \\    --nx N --ny N
        \\    --width L --height L
        \\    --courant C
        \\    --dt DT
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

fn print3DUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli euler reference [options]
        \\
        \\  options:
        \\    --steps N
        \\    --n0 N --n1 N --n2 N
        \\    --extent0 L --extent1 L --extent2 L
        \\    --courant C
        \\    --dt DT
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
