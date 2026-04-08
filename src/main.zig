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
    const stderr = (std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter();
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
