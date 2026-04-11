const std = @import("std");
const commands = @import("example_commands");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const argv = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, argv);
    try run(allocator, argv);
}

pub fn run(allocator: std.mem.Allocator, argv: []const [:0]const u8) !void {
    if (argv.len < 2) {
        printUsage();
        std.process.exit(2);
    }

    const cmd = argv[1];
    if (eql(cmd, "list")) {
        listCommands();
        return;
    }
    if (eql(cmd, "--help") or eql(cmd, "-h") or eql(cmd, "help")) {
        printUsage();
        return;
    }

    inline for (commands.subcommands) |sub| {
        if (eql(cmd, sub.name)) {
            return sub.run(allocator, argv[1..]);
        }
    }

    std.debug.print("error: unknown subcommand '{s}'\n\n", .{cmd});
    printUsage();
    std.process.exit(2);
}

pub fn printUsage() void {
    std.debug.print(
        \\
        \\  flux — discrete exterior calculus simulation suite
        \\
        \\  usage: flux <subcommand> [options]
        \\         flux list
        \\         flux <subcommand> --help
        \\
        \\  subcommands:
        \\
    , .{});
    inline for (commands.subcommands) |sub| {
        std.debug.print("    {s:<20}  {s}\n", .{ sub.name, sub.summary });
    }
    std.debug.print(
        \\
        \\  ask a family for its full flag list:
        \\    flux maxwell --dim 2 --help
        \\    flux euler --dim 3 --help
        \\    flux diffusion --help
        \\
        \\  examples:
        \\    zig build run -- diffusion --surface plane --grid 32 --frames 4
        \\    zig build run -- maxwell --dim 2 --demo cavity --frames 50
        \\    zig build run -- euler --dim 3 --steps 100 --frames 5
        \\
    , .{});
}

fn listCommands() void {
    inline for (commands.subcommands) |sub| {
        std.debug.print("{s}\n", .{sub.name});
    }
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
