const std = @import("std");
const common = @import("examples_common");
const maxwell_2d = @import("maxwell_2d");
const maxwell_3d = @import("maxwell_3d");
const euler_2d = @import("euler_2d");
const euler_3d = @import("euler_3d");
const diffusion = @import("diffusion");

const DiffusionSurface = enum { plane, sphere };
const UsageFn = *const fn () void;

pub const Subcommand = struct {
    name: []const u8,
    summary: []const u8,
    run: *const fn (allocator: std.mem.Allocator, args: []const [:0]const u8) anyerror!void,
};

pub const subcommands = [_]Subcommand{
    .{ .name = "maxwell-2d", .summary = "2D Maxwell on simplicial meshes (cavity, dipole)", .run = runMaxwell2d },
    .{ .name = "maxwell-3d", .summary = "3D Maxwell TM_110 cavity on tetrahedra", .run = runMaxwell3d },
    .{ .name = "euler-2d", .summary = "2D incompressible Euler vorticity-stream", .run = runEuler2d },
    .{ .name = "euler-3d", .summary = "3D incompressible Euler with helicity conservation", .run = runEuler3d },
    .{ .name = "diffusion", .summary = "Scalar diffusion on a plane or curved surface", .run = runDiffusion },
};

const usage_maxwell_2d =
    \\  flux-examples maxwell-2d — 2D electromagnetic simulation on simplicial meshes
    \\
    \\  usage: flux-examples maxwell-2d [--demo <name>] [options]
    \\
    \\  demos:
    \\    dipole    (default) point dipole radiating in a PEC cavity
    \\    cavity    TE10 standing wave — exact mode, source-free
    \\
    \\  mesh & physics:
    \\    --grid N          cells per side (default: 32)
    \\    --domain L        domain side length (default: 1.0)
    \\    --courant C       Courant number, dt = C·h (default: 0.1)
    \\    --dt DT           explicit timestep override (computes courant = dt/h)
    \\
    \\  dipole source (ignored for cavity):
    \\    --frequency F     source frequency in Hz (default: TE10 resonance)
    \\    --amplitude A     source amplitude (default: 1.0)
    \\
    \\  time stepping:
    \\    --steps N         number of timesteps (default: 1000)
    \\
    \\  output:
    \\    --output DIR      VTK output directory (default: output)
    \\    --frames N        number of snapshots (default: 100)
    \\
;

const usage_maxwell_3d =
    \\  flux-examples maxwell-3d — TM_110 cavity resonance on tetrahedra
    \\
    \\  usage: flux-examples maxwell-3d [options]
    \\
    \\  options:
    \\    --nx N                tetrahedral cells in x (default: 2)
    \\    --ny N                tetrahedral cells in y (default: 2)
    \\    --nz N                tetrahedral cells in z (default: 2)
    \\    --width L             cavity width  (default: 1.0)
    \\    --height L            cavity height (default: 1.0)
    \\    --depth L             cavity depth  (default: 1.0)
    \\    --steps N             leapfrog steps (default: 1000)
    \\    --dt DT               fixed timestep (default: 0.01)
    \\    --output DIR          write VTK snapshots into DIR
    \\    --output-interval N   write every N steps (overrides --frames)
    \\    --frames N            number of evenly-spaced snapshots
    \\
;

const usage_euler_2d =
    \\  flux-examples euler-2d — 2D incompressible Euler vorticity-stream
    \\
    \\  usage: flux-examples euler-2d [options]
    \\
    \\  options:
    \\    --demo NAME       gaussian (invariant reference) or dipole (dynamic showcase)
    \\    --grid N          grid cells per side (default: 16)
    \\    --steps N         number of timesteps (default: 1000)
    \\    --domain L        square domain side length (default: 1.0)
    \\    --cfl C           timestep scale dt = C·h (default: 0.1)
    \\    --dt DT           explicit timestep override (computes cfl = dt/h)
    \\    --output DIR      output directory for VTK files (default: output/euler_2d)
    \\    --frames N        number of snapshots to write (default: 50)
    \\
;

const usage_euler_3d =
    \\  flux-examples euler-3d — 3D incompressible Euler with helicity
    \\
    \\  usage: flux-examples euler-3d [options]
    \\
    \\  options:
    \\    --steps N             number of timesteps (default: 1000)
    \\    --nx N                x-axis cube count (default: 2)
    \\    --ny N                y-axis cube count (default: 2)
    \\    --nz N                z-axis cube count (default: 2)
    \\    --width L             domain width (default: 1.0)
    \\    --height L            domain height (default: 1.0)
    \\    --depth L             domain depth (default: 1.0)
    \\    --dt DT               timestep size (default: 0.01)
    \\    --output DIR          output directory for VTK snapshots
    \\    --output-interval N   snapshot cadence in steps (overrides --frames)
    \\    --frames N            number of evenly-spaced snapshots
    \\
;

const usage_diffusion =
    \\  flux-examples diffusion — scalar diffusion on a plane or sphere
    \\
    \\  usage: flux-examples diffusion [options]
    \\
    \\  options:
    \\    --surface NAME     plane (default) or sphere
    \\    --steps N          backward-Euler steps
    \\    --dt-scale C       diffusion timestep scale hint
    \\    --dt DT            explicit timestep override
    \\    --output DIR       output directory for VTK files
    \\    --frames N         number of snapshots to write
    \\
    \\  plane-only:
    \\    --grid N           grid cells per side (default: 8)
    \\    --domain L         square domain side length (default: 1.0)
    \\
    \\  sphere-only:
    \\    --refinement N     spherical mesh refinement level (default: 0)
    \\    --final-time T     total simulated time (default: 0.05)
    \\
;

pub fn runMaxwell2d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = maxwell_2d.Config{};
    var shared = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (eql(arg, "--demo")) {
            const value = try parser.requireValue("--demo");
            if (eql(value, "dipole")) {
                config.demo = .dipole;
            } else if (eql(value, "cavity")) {
                config.demo = .cavity;
            } else {
                std.debug.print("error: unknown demo '{s}'. available: dipole, cavity\n", .{value});
                std.process.exit(2);
            }
            continue;
        }
        if (eql(arg, "--courant")) {
            config.courant = try parser.parseF64("--courant");
            continue;
        }
        if (eql(arg, "--frequency")) {
            config.frequency = try parser.parseF64("--frequency");
            continue;
        }
        if (eql(arg, "--amplitude")) {
            config.amplitude = try parser.parseF64("--amplitude");
            continue;
        }
        if (try parser.tryCommon(arg, &shared)) continue;
        failUnknownFlag(arg, printMaxwell2dUsage);
    }

    if (shared.help) {
        printMaxwell2dUsage();
        return;
    }

    common.applySharedFields(&config, shared);
    if (shared.dt) |dt_value| {
        const h = config.domain / @as(f64, @floatFromInt(config.grid));
        config.courant = dt_value / h;
    }

    try maxwell_2d.runDriver(allocator, config);
}

pub fn runMaxwell3d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    const config = try parseBox3Command(maxwell_3d.Config, args, printMaxwell3dUsage, applyMaxwell3dShared) orelse return;
    try maxwell_3d.runDriver(allocator, config);
}

pub fn runEuler2d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = euler_2d.Config{};
    var shared = common.Common{};
    var output_explicit = false;
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (eql(arg, "--demo")) {
            const value = try parser.requireValue("--demo");
            if (eql(value, "gaussian")) {
                config.demo = .gaussian;
            } else if (eql(value, "dipole")) {
                config.demo = .dipole;
            } else {
                std.debug.print("error: unknown demo '{s}'. available: gaussian, dipole\n", .{value});
                std.process.exit(2);
            }
            continue;
        }
        if (eql(arg, "--cfl")) {
            config.cfl = try parser.parseF64("--cfl");
            continue;
        }
        if (try parser.tryCommon(arg, &shared)) {
            if (shared.output_dir != null) output_explicit = true;
            continue;
        }
        failUnknownFlag(arg, printEuler2dUsage);
    }

    if (shared.help) {
        printEuler2dUsage();
        return;
    }

    common.applySharedFields(&config, shared);
    if (shared.dt) |dt_value| {
        const h = config.domain / @as(f64, @floatFromInt(config.grid));
        config.cfl = dt_value / h;
    }
    if (!output_explicit and config.demo == .dipole) {
        config.output_dir = "output/euler_dipole";
    }

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try euler_2d.run(allocator, config, stderr);
    const drift = @abs(result.circulation_final - result.circulation_initial);
    try stderr.print(
        "elapsed={d:.3}s snapshots={d} drift={e} output={s}\n",
        .{ result.elapsed_s, result.snapshot_count, drift, config.output_dir },
    );
}

pub fn runEuler3d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    const config = try parseBox3Command(euler_3d.Config, args, printEuler3dUsage, applyEuler3dShared) orelse return;

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try euler_3d.run(allocator, config, stderr);
    try stderr.print(
        "elapsed={d:.3}s snapshots={d} drift={e}\n",
        .{
            result.elapsed_s,
            result.snapshot_count,
            @abs(result.helicity_final - result.helicity_initial),
        },
    );
}

pub fn runDiffusion(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var surface_kind: DiffusionSurface = .plane;
    var plane_config = diffusion.PlaneConfig{};
    var sphere_config = diffusion.SphereConfig{};
    var shared = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (eql(arg, "--surface")) {
            const value = try parser.requireValue("--surface");
            if (eql(value, "plane")) {
                surface_kind = .plane;
            } else if (eql(value, "sphere")) {
                surface_kind = .sphere;
            } else {
                std.debug.print("error: unknown surface '{s}'. available: plane, sphere\n", .{value});
                std.process.exit(2);
            }
            continue;
        }
        if (try parser.tryCommon(arg, &shared)) continue;
        if (eql(arg, "--dt-scale")) {
            const value = try parser.parseF64("--dt-scale");
            plane_config.dt_scale = value;
            sphere_config.dt_scale = value;
            continue;
        }
        failUnknownFlag(arg, printDiffusionUsage);
    }

    if (shared.help) {
        printDiffusionUsage();
        return;
    }

    switch (surface_kind) {
        .plane => {
            common.applySharedFields(&plane_config, shared);
            if (shared.dt) |value| plane_config.dt_override = value;
            try runPlaneDiffusion(allocator, plane_config);
        },
        .sphere => {
            common.applySharedFields(&sphere_config, shared);
            if (shared.dt) |dt_value| {
                sphere_config.final_time = dt_value * @as(f64, @floatFromInt(sphere_config.steps));
            }
            try runSphereDiffusion(allocator, sphere_config);
        },
    }
}

fn runPlaneDiffusion(allocator: std.mem.Allocator, config: diffusion.PlaneConfig) !void {
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try diffusion.runPlane(allocator, config, stderr);
    try stderr.print(
        "elapsed={d:.3}s steps={d} snapshots={d} l2_error={e} output={s}\n",
        .{ result.elapsed_s, result.steps, result.snapshot_count, result.l2_error, config.output_dir },
    );
}

fn runSphereDiffusion(allocator: std.mem.Allocator, config: diffusion.SphereConfig) !void {
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    _ = try diffusion.runSphere(allocator, config, stderr);
}

fn parseBox3Command(
    comptime ConfigType: type,
    args: []const [:0]const u8,
    usage_fn: UsageFn,
    apply_shared: *const fn (*ConfigType, common.Common) void,
) !?ConfigType {
    var config = ConfigType{};
    var shared = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (try common.tryBox3Flag(&parser, arg, &config)) continue;
        if (try parser.tryCommon(arg, &shared)) continue;
        failUnknownFlag(arg, usage_fn);
    }

    if (shared.help) {
        usage_fn();
        return null;
    }

    common.applySharedFields(&config, shared);
    apply_shared(&config, shared);
    return config;
}

fn applyFixedDtWithFrameInterval(cfg: anytype, shared: common.Common) void {
    if (shared.dt) |value| cfg.dt = value;
    if (shared.frames) |frames| {
        cfg.output_interval = common.framesToInterval(cfg.steps, frames);
    }
}

fn applyMaxwell3dShared(cfg: *maxwell_3d.Config, shared: common.Common) void {
    applyFixedDtWithFrameInterval(cfg, shared);
}

fn applyEuler3dShared(cfg: *euler_3d.Config, shared: common.Common) void {
    applyFixedDtWithFrameInterval(cfg, shared);
}

fn failUnknownFlag(arg: []const u8, usage_fn: UsageFn) noreturn {
    std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
    usage_fn();
    std.process.exit(2);
}

fn printUsageText(text: []const u8) void {
    std.debug.print("{s}", .{text});
}

fn printMaxwell2dUsage() void {
    printUsageText(usage_maxwell_2d);
}

fn printMaxwell3dUsage() void {
    printUsageText(usage_maxwell_3d);
}

fn printEuler2dUsage() void {
    printUsageText(usage_euler_2d);
}

fn printEuler3dUsage() void {
    printUsageText(usage_euler_3d);
}

fn printDiffusionUsage() void {
    printUsageText(usage_diffusion);
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
