const std = @import("std");
const diffusion = @import("diffusion.zig");
const euler = @import("euler.zig");
const maxwell = @import("maxwell.zig");

pub const Subcommand = struct {
    name: []const u8,
    summary: []const u8,
    run: *const fn (allocator: std.mem.Allocator, args: []const [:0]const u8) anyerror!void,
};

pub const subcommands = [_]Subcommand{
    .{ .name = diffusion.name, .summary = diffusion.summary, .run = diffusion.run },
    .{ .name = euler.name, .summary = euler.summary, .run = euler.run },
    .{ .name = maxwell.name, .summary = maxwell.summary, .run = maxwell.run },
};

test {
    std.testing.refAllDeclsRecursive(@This());
}
