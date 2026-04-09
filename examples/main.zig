//! Single entrypoint for the flux example suite.
//!
//! Every demo lives behind a subcommand of the `flux-examples` binary so
//! users only need to learn one CLI surface:
//!
//!     flux-examples list
//!     flux-examples maxwell-2d --steps 1000 --frames 50
//!     flux-examples heat --grid 32 --frames 4
//!     flux-examples euler-3d --nx 4 --ny 4 --nz 4
//!
//! Each subcommand parses its own physics-specific flags on top of the shared
//! `examples_common.cli.Common` schema (`--steps`, `--dt`, `--output`,
//! `--frames`, `--grid`, `--domain`, `--refinement`, `--final-time`).

const std = @import("std");
const registry = @import("examples_registry");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const argv = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, argv);

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

    inline for (registry.subcommands) |sub| {
        if (eql(cmd, sub.name)) {
            return sub.run(allocator, argv[1..]);
        }
    }

    std.debug.print("error: unknown subcommand '{s}'\n\n", .{cmd});
    printUsage();
    std.process.exit(2);
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples — discrete exterior calculus simulation suite
        \\
        \\  usage: flux-examples <subcommand> [options]
        \\         flux-examples list
        \\         flux-examples <subcommand> --help
        \\
        \\  subcommands:
        \\
    , .{});
    inline for (registry.subcommands) |sub| {
        std.debug.print("    {s:<20}  {s}\n", .{ sub.name, sub.summary });
    }
    std.debug.print(
        \\
        \\  shared flags (recognized by every subcommand):
        \\    --steps N         number of timesteps
        \\    --dt DT           override the physics-derived timestep
        \\    --output DIR      output directory for VTK snapshots
        \\    --frames N        number of evenly-spaced snapshots to write
        \\    --grid N          structured grid cells per side (where applicable)
        \\    --domain L        domain side length (where applicable)
        \\    --refinement N    refinement level (where applicable)
        \\    --final-time T    final simulated time (where applicable)
        \\    --help, -h        show subcommand-specific help
        \\
        \\  build:
        \\    zig build examples           # build the umbrella binary
        \\    zig build run-maxwell-2d -- --frames 50
        \\    zig build run-heat -- --grid 32
        \\
    , .{});
}

fn listCommands() void {
    inline for (registry.subcommands) |sub| {
        std.debug.print("{s}\n", .{sub.name});
    }
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
