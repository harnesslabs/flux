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

const maxwell_2d_cli = @import("maxwell_2d/cli.zig");
const maxwell_3d_cli = @import("maxwell_3d/cli.zig");
const euler_2d_cli = @import("euler_2d/cli.zig");
const euler_3d_cli = @import("euler_3d/cli.zig");
const heat_cli = @import("heat/cli.zig");
const diffusion_surface_cli = @import("diffusion_surface/cli.zig");

/// One row of the dispatch table. The `run` function receives the entire
/// argv tail (subcommand name plus its own flags) and is responsible for
/// printing its own usage on `--help`.
const Subcommand = struct {
    name: []const u8,
    summary: []const u8,
    run: *const fn (allocator: std.mem.Allocator, args: []const [:0]const u8) anyerror!void,
};

const subcommands = [_]Subcommand{
    .{
        .name = "maxwell-2d",
        .summary = "2D Maxwell on simplicial meshes (cavity, dipole)",
        .run = maxwell_2d_cli.run,
    },
    .{
        .name = "maxwell-3d",
        .summary = "3D Maxwell TM_110 cavity on tetrahedra",
        .run = maxwell_3d_cli.run,
    },
    .{
        .name = "euler-2d",
        .summary = "2D incompressible Euler vorticity-stream",
        .run = euler_2d_cli.run,
    },
    .{
        .name = "euler-3d",
        .summary = "3D incompressible Euler with helicity conservation",
        .run = euler_3d_cli.run,
    },
    .{
        .name = "heat",
        .summary = "Implicit heat equation via backward Euler + CG",
        .run = heat_cli.run,
    },
    .{
        .name = "diffusion-surface",
        .summary = "Heat equation on a curved surface (Riemannian Hodge)",
        .run = diffusion_surface_cli.run,
    },
};

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

    inline for (subcommands) |sub| {
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
    inline for (subcommands) |sub| {
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
    inline for (subcommands) |sub| {
        std.debug.print("{s}\n", .{sub.name});
    }
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
