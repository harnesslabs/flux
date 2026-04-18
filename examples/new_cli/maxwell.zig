const std = @import("std");
const parse = @import("parse.zig");
const maxwell = @import("new_maxwell");

pub const name = "maxwell";
pub const summary = "New Maxwell example CLI over examples/new_maxwell";

const Scenario = enum {
    dipole,
    cavity,
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

    const scenario: Scenario = if (parse.eql(scenario_arg, "dipole"))
        .dipole
    else if (parse.eql(scenario_arg, "cavity"))
        .cavity
    else {
        std.debug.print("error: unknown Maxwell scenario '{s}'\n\n", .{scenario_arg});
        printUsage();
        std.process.exit(2);
    };

    switch (scenario) {
        .dipole => try runDipole(allocator, &parser),
        .cavity => try runCavity(allocator, &parser),
    }
}

fn runDipole(allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = maxwell.DipoleConfig2D{};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            printDipoleUsage();
            return;
        }
        if (parse.eql(arg, "--steps")) {
            config.steps = try parser.parseU32("--steps");
            continue;
        }
        if (parse.eql(arg, "--nx")) {
            const value = try parser.parseU32("--nx");
            config.counts[0] = value;
            continue;
        }
        if (parse.eql(arg, "--ny")) {
            const value = try parser.parseU32("--ny");
            config.counts[1] = value;
            continue;
        }
        if (parse.eql(arg, "--width")) {
            const value = try parser.parseF64("--width");
            config.extents[0] = value;
            continue;
        }
        if (parse.eql(arg, "--height")) {
            const value = try parser.parseF64("--height");
            config.extents[1] = value;
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
        if (parse.eql(arg, "--frequency")) {
            config.frequency_hz = try parser.parseF64("--frequency");
            continue;
        }
        if (parse.eql(arg, "--amplitude")) {
            config.amplitude = try parser.parseF64("--amplitude");
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
        failUnknownFlag(arg, printDipoleUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = try maxwell.runDipole(allocator, config, stderr);
    try stderr.print(
        "dipole: elapsed={d:.3}s snapshots={d} energy={e}\n",
        .{ result.elapsed_s, result.snapshot_count, result.summary.energy_final },
    );
}

fn runCavity(allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var dim: u8 = 2;
    if (parser.peek()) |arg| {
        if (parse.eql(arg, "--dim")) {
            _ = parser.next();
            dim = @intCast(try parser.parseU32("--dim"));
        }
    }

    switch (dim) {
        2 => try runCavityDim(2, allocator, parser),
        3 => try runCavityDim(3, allocator, parser),
        else => {
            std.debug.print("error: Maxwell cavity only supports --dim 2 or --dim 3\n", .{});
            std.process.exit(2);
        },
    }
}

fn runCavityDim(comptime dim: u8, allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = maxwell.CavityConfig(dim){};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            printCavityUsage();
            return;
        }
        if (parse.eql(arg, "--steps")) {
            config.steps = try parser.parseU32("--steps");
            continue;
        }
        if (parse.eql(arg, "--reference")) {
            config.reference = true;
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
        if (try consumeCounts(dim, parser, arg, &config.counts)) continue;
        if (try consumeExtents(dim, parser, arg, &config.extents)) continue;
        failUnknownFlag(arg, printCavityUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = try maxwell.runCavity(dim, allocator, config, stderr);
    if (config.reference) {
        try stderr.print(
            "cavity_{d}d: elapsed={d:.3}s snapshots={d} energy={e} electric_l2={e} magnetic_l2={e}\n",
            .{
                dim,
                result.elapsed_s,
                result.snapshot_count,
                result.summary.energy_final,
                result.summary.electric_l2_final.?,
                result.summary.magnetic_l2_final.?,
            },
        );
    } else {
        try stderr.print(
            "cavity_{d}d: elapsed={d:.3}s snapshots={d} energy={e}\n",
            .{ dim, result.elapsed_s, result.snapshot_count, result.summary.energy_final },
        );
    }
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
        \\  flux-new-cli maxwell <scenario> [options]
        \\
        \\  scenarios:
        \\    dipole    2D point dipole radiation
        \\    cavity    2D or 3D cavity mode
        \\
        \\  ask a scenario for details:
        \\    zig build run-new-cli -- maxwell dipole --help
        \\    zig build run-new-cli -- maxwell cavity --help
        \\
    , .{});
}

fn printDipoleUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli maxwell dipole [options]
        \\
        \\  options:
        \\    --steps N
        \\    --nx N --ny N
        \\    --width L --height L
        \\    --courant C
        \\    --dt DT
        \\    --frequency F
        \\    --amplitude A
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

fn printCavityUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli maxwell cavity [--dim 2|3] [options]
        \\
        \\  options:
        \\    --steps N
        \\    --n0 N --n1 N [--n2 N]
        \\    --extent0 L --extent1 L [--extent2 L]
        \\    --courant C
        \\    --dt DT
        \\    --reference
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
