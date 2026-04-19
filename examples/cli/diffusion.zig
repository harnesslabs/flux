const std = @import("std");
const diffusion = @import("diffusion");
const parse = @import("parse.zig");

pub const name = "diffusion";
pub const summary = "Diffusion examples";

const Scenario = enum {
    plane,
    sphere,
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

    const scenario: Scenario = if (parse.eql(scenario_arg, "plane"))
        .plane
    else if (parse.eql(scenario_arg, "sphere"))
        .sphere
    else {
        std.debug.print("error: unknown diffusion scenario '{s}'\n\n", .{scenario_arg});
        printUsage();
        std.process.exit(2);
    };

    switch (scenario) {
        .plane => try runPlane(allocator, &parser),
        .sphere => try runSphere(allocator, &parser),
    }
}

fn runPlane(allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = diffusion.PlaneConfig{};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            printPlaneUsage();
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
        if (parse.eql(arg, "--dt-scale")) {
            config.dt_scale = try parser.parseF64("--dt-scale");
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
        failUnknownFlag(arg, printPlaneUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = try diffusion.runPlane(allocator, config, stderr);
    try stderr.print(
        "diffusion_plane: elapsed={d:.3}s snapshots={d} l2_error={e}\n",
        .{ result.elapsed_s, result.snapshot_count, result.summary.l2_error_final },
    );
}

fn runSphere(allocator: std.mem.Allocator, parser: *parse.Parser) !void {
    var config = diffusion.SphereConfig{};

    while (parser.next()) |arg| {
        if (parse.eql(arg, "--help") or parse.eql(arg, "-h")) {
            printSphereUsage();
            return;
        }
        if (parse.eql(arg, "--steps")) {
            config.steps = try parser.parseU32("--steps");
            continue;
        }
        if (parse.eql(arg, "--refinement")) {
            config.refinement = try parser.parseU32("--refinement");
            continue;
        }
        if (parse.eql(arg, "--radius")) {
            config.radius = try parser.parseF64("--radius");
            continue;
        }
        if (parse.eql(arg, "--final-time")) {
            config.final_time = try parser.parseF64("--final-time");
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
        failUnknownFlag(arg, printSphereUsage);
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    const result = try diffusion.runSphere(allocator, config, stderr);
    try stderr.print(
        "diffusion_sphere: elapsed={d:.3}s snapshots={d} l2_error={e}\n",
        .{ result.elapsed_s, result.snapshot_count, result.summary.l2_error_final },
    );
}

fn failUnknownFlag(arg: []const u8, usage: *const fn () void) noreturn {
    std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
    usage();
    std.process.exit(2);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli diffusion <scenario> [options]
        \\
        \\  scenarios:
        \\    plane     planar diffusion with exact sine-mode reference
        \\    sphere    surface diffusion on a sphere
        \\
        \\  ask a scenario for details:
        \\    zig build run-new-cli -- diffusion plane --help
        \\    zig build run-new-cli -- diffusion sphere --help
        \\
    , .{});
}

fn printPlaneUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli diffusion plane [options]
        \\
        \\  options:
        \\    --steps N
        \\    --nx N --ny N
        \\    --width L --height L
        \\    --dt-scale C
        \\    --dt DT
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

fn printSphereUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli diffusion sphere [options]
        \\
        \\  options:
        \\    --steps N
        \\    --refinement N
        \\    --radius R
        \\    --final-time T
        \\    --output DIR
        \\    --frames N
        \\    --output-interval N
        \\
    , .{});
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
