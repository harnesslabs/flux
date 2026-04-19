const std = @import("std");
const commands = @import("cli_commands");

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

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-new-cli <subcommand> [options]
        \\
        \\  subcommands:
        \\
    , .{});
    inline for (commands.subcommands) |sub| {
        std.debug.print("    {s:<16}  {s}\n", .{ sub.name, sub.summary });
    }
    std.debug.print("\n", .{});
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
