const std = @import("std");
const heat = @import("heat.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    const config = parseArgs(args) catch return;
    const stderr = (std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter();
    const result = try heat.run(allocator, config, stderr);
    try stderr.print(
        "elapsed={d:.3}s steps={d} snapshots={d} l2_error={e} output={s}\n",
        .{ result.elapsed_s, result.steps, result.snapshot_count, result.l2_error, config.output_dir },
    );
}

test {
    _ = @import("heat.zig");
}

const ParseError = error{InvalidArgument};

fn parseArgs(args: []const [:0]const u8) ParseError!heat.Config {
    var config = heat.Config{};
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (eql(arg, "--help") or eql(arg, "-h")) {
            printUsage();
            std.process.exit(0);
        } else if (eql(arg, "--grid")) {
            config.grid = parseU32(args, &i, "--grid") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--steps")) {
            config.steps = parseU32(args, &i, "--steps") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--domain")) {
            config.domain = parseF64(args, &i, "--domain") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--output")) {
            config.output_dir = nextArg(args, &i) orelse return flagError("--output");
        } else if (eql(arg, "--frames")) {
            config.frames = parseU32(args, &i, "--frames") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--dt-scale")) {
            config.dt_scale = parseF64(args, &i, "--dt-scale") orelse return ParseError.InvalidArgument;
        } else {
            std.debug.print("error: unknown argument '{s}'\n\n", .{arg});
            printUsage();
            return ParseError.InvalidArgument;
        }
    }
    return config;
}

fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn nextArg(args: []const [:0]const u8, i: *usize) ?[]const u8 {
    if (i.* + 1 < args.len) {
        i.* += 1;
        return args[i.*];
    }
    return null;
}

fn flagError(flag: []const u8) ParseError {
    std.debug.print("error: {s} requires a value\n", .{flag});
    return ParseError.InvalidArgument;
}

fn parseU32(args: []const [:0]const u8, i: *usize, flag: []const u8) ?u32 {
    const value = nextArg(args, i) orelse {
        std.debug.print("error: {s} requires a value\n", .{flag});
        return null;
    };
    return std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("error: invalid {s} value: {s}\n", .{ flag, value });
        return null;
    };
}

fn parseF64(args: []const [:0]const u8, i: *usize, flag: []const u8) ?f64 {
    const value = nextArg(args, i) orelse {
        std.debug.print("error: {s} requires a value\n", .{flag});
        return null;
    };
    return std.fmt.parseFloat(f64, value) catch {
        std.debug.print("error: invalid {s} value: {s}\n", .{ flag, value });
        return null;
    };
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux — implicit heat equation example
        \\
        \\  usage: zig build -Doptimize=ReleaseFast example-heat -- [options]
        \\
        \\  options:
        \\    --grid N          grid cells per side (default: 8)
        \\    --steps N         backward-Euler steps (default: 8)
        \\    --domain L        square domain side length (default: 1.0)
        \\    --dt-scale C      timestep scale dt = C·h² (default: 0.1)
        \\    --output DIR      output directory for VTK files (default: output/heat)
        \\    --frames N        number of snapshots to write (default: 4)
        \\
    , .{});
}
