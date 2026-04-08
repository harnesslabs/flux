//! flux CLI entry point.
//!
//! The flux library provides composable DEC operators for PDE simulation.
//! Physics-specific simulations live in `examples/` and consume flux as
//! a library dependency.
//!
//! Available examples:
//!   zig build run-maxwell-2d -- [options]   # 2D electromagnetic simulation

const std = @import("std");

pub fn main() void {
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    stderr.writeAll(
        \\
        \\  flux — composable, type-safe PDE solver framework
        \\
        \\  This binary is a placeholder. Physics simulations live in examples/:
        \\
        \\    zig build run-maxwell-2d -- [options]   # 2D Maxwell simulation
        \\
        \\  Run `zig build --help` to see all available build steps.
        \\
        \\
    ) catch {};
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
