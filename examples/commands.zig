const std = @import("std");
const common = @import("examples_common");
const maxwell = @import("maxwell");
const euler = @import("euler");
const diffusion = @import("diffusion");

const UsageFn = *const fn () void;

pub const Subcommand = struct {
    name: []const u8,
    summary: []const u8,
    run: *const fn (allocator: std.mem.Allocator, args: []const [:0]const u8) anyerror!void,
};

pub const subcommands = [_]Subcommand{
    .{ .name = "maxwell", .summary = "Discrete Maxwell simulation in 2D or 3D", .run = runMaxwell },
    .{ .name = "euler", .summary = "Incompressible Euler simulation in 2D or 3D", .run = runEuler },
    .{ .name = "diffusion", .summary = "Scalar diffusion on a plane or sphere", .run = runDiffusion },
};

const usage_maxwell_2d =
    \\  flux maxwell --dim 2 [--demo <name>] [options]
    \\
    \\  demos:
    \\    dipole    (default) point dipole radiating in a PEC cavity
    \\    cavity    TE10 standing wave, source-free
    \\
    \\  options:
    \\    --grid N          cells per side (default: 32)
    \\    --domain L        domain side length (default: 1.0)
    \\    --courant C       Courant number, dt = C*h (default: 0.1)
    \\    --dt DT           explicit timestep override
    \\    --frequency F     source frequency in Hz
    \\    --amplitude A     source amplitude
    \\    --steps N         number of timesteps
    \\    --output DIR      VTK output directory
    \\    --frames N        number of snapshots
    \\
;

const usage_maxwell_3d =
    \\  flux maxwell --dim 3 [options]
    \\
    \\  options:
    \\    --nx N --ny N --nz N   tetrahedral cells per axis
    \\    --width L              cavity width
    \\    --height L             cavity height
    \\    --depth L              cavity depth
    \\    --steps N              leapfrog steps
    \\    --dt DT                fixed timestep
    \\    --output DIR           write VTK snapshots into DIR
    \\    --output-interval N    write every N steps
    \\    --frames N             evenly-spaced snapshots
    \\
;

const usage_euler_2d =
    \\  flux euler --dim 2 [options]
    \\
    \\  options:
    \\    --demo NAME       gaussian (default) or dipole
    \\    --grid N          grid cells per side (default: 16)
    \\    --steps N         number of timesteps (default: 1000)
    \\    --domain L        square domain side length (default: 1.0)
    \\    --cfl C           timestep scale dt = C*h (default: 0.1)
    \\    --dt DT           explicit timestep override
    \\    --output DIR      output directory for VTK files
    \\    --frames N        number of snapshots
    \\
;

const usage_euler_3d =
    \\  flux euler --dim 3 [options]
    \\
    \\  options:
    \\    --steps N             number of timesteps
    \\    --nx N --ny N --nz N  cube count per axis
    \\    --width L             domain width
    \\    --height L            domain height
    \\    --depth L             domain depth
    \\    --dt DT               timestep size
    \\    --output DIR          output directory for VTK snapshots
    \\    --output-interval N   snapshot cadence in steps
    \\    --frames N            evenly-spaced snapshots
    \\
;

const usage_diffusion =
    \\  flux diffusion [options]
    \\
    \\  options:
    \\    --surface NAME     plane (default) or sphere
    \\    --steps N          backward-Euler steps
    \\    --dt-scale C       diffusion timestep scale hint
    \\    --dt DT            explicit timestep override
    \\    --output DIR       output directory for VTK files
    \\    --frames N         number of snapshots to write
    \\
    \\  plane:
    \\    --grid N           grid cells per side
    \\    --domain L         square domain side length
    \\
    \\  sphere:
    \\    --refinement N     spherical mesh refinement level
    \\    --final-time T     total simulated time
    \\
;

pub fn runMaxwell(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    try runByDimension(allocator, args, runMaxwell2d, runMaxwell3d);
}

pub fn runEuler(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    try runByDimension(allocator, args, runEuler2d, runEuler3d);
}

fn runMaxwell2d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = maxwell.Config(2){};
    var shared = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (try consumeDimensionFlag(&parser, arg)) continue;
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

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    _ = try maxwell.run(2, allocator, config, stderr);
}

fn runMaxwell3d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    const config = try parseBox3Command(
        maxwell.Config(3),
        args,
        printMaxwell3dUsage,
        applyFixedDtWithFrameInterval,
    ) orelse return;
    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    _ = try maxwell.run(3, allocator, config, stderr);
}

fn runEuler2d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = euler.Config(2){};
    var shared = common.Common{};
    var output_explicit = false;
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (try consumeDimensionFlag(&parser, arg)) continue;
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
    const result = try euler.run(2, allocator, config, stderr);
    const drift = @abs(result.circulation_final - result.circulation_initial);
    try stderr.print(
        "elapsed={d:.3}s snapshots={d} drift={e} output={s}\n",
        .{ result.elapsed_s, result.snapshot_count, drift, config.output_dir },
    );
}

fn runEuler3d(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    const config = try parseBox3Command(
        euler.Config(3),
        args,
        printEuler3dUsage,
        applyFixedDtWithFrameInterval,
    ) orelse return;

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;
    const result = try euler.run(3, allocator, config, stderr);
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
    var surface_kind: diffusion.SurfaceKind = .plane;
    var plane_config = diffusion.Config(.plane){};
    var sphere_config = diffusion.Config(.sphere){};
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
            var stderr_buffer: [1024]u8 = undefined;
            var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
            const stderr = &stderr_writer.interface;
            const result = try diffusion.run(.plane, allocator, plane_config, stderr);
            try stderr.print(
                "elapsed={d:.3}s steps={d} snapshots={d} l2_error={e} output={s}\n",
                .{ result.elapsed_s, result.steps, result.snapshot_count, result.l2_error, plane_config.output_dir },
            );
        },
        .sphere => {
            common.applySharedFields(&sphere_config, shared);
            if (shared.dt) |dt_value| {
                sphere_config.final_time = dt_value * @as(f64, @floatFromInt(sphere_config.steps));
            }
            var stderr_buffer: [1024]u8 = undefined;
            var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
            const stderr = &stderr_writer.interface;
            _ = try diffusion.run(.sphere, allocator, sphere_config, stderr);
        },
    }
}

fn runByDimension(
    allocator: std.mem.Allocator,
    args: []const [:0]const u8,
    run2: *const fn (std.mem.Allocator, []const [:0]const u8) anyerror!void,
    run3: *const fn (std.mem.Allocator, []const [:0]const u8) anyerror!void,
) !void {
    return switch (try parseDimensionArg(args, 2)) {
        2 => try run2(allocator, args),
        3 => try run3(allocator, args),
        else => unreachable,
    };
}

fn parseDimensionArg(args: []const [:0]const u8, default_dimension: u8) !u8 {
    var parser = common.Parser.init(args);
    _ = parser.next();
    while (parser.next()) |arg| {
        if (!eql(arg, "--dim")) continue;
        const value = try parser.parseU32("--dim");
        return switch (value) {
            2, 3 => @intCast(value),
            else => {
                std.debug.print("error: unsupported --dim {d}; expected 2 or 3\n", .{value});
                std.process.exit(2);
            },
        };
    }
    return default_dimension;
}

fn parseBox3Command(
    comptime ConfigType: type,
    args: []const [:0]const u8,
    usage_fn: UsageFn,
    comptime apply_shared: anytype,
) !?ConfigType {
    var config = ConfigType{};
    var shared = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (try consumeDimensionFlag(&parser, arg)) continue;
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

fn consumeDimensionFlag(parser: *common.Parser, arg: []const u8) !bool {
    if (!eql(arg, "--dim")) return false;
    _ = try parser.requireValue("--dim");
    return true;
}

fn applyFixedDtWithFrameInterval(cfg: anytype, shared: common.Common) void {
    if (shared.dt) |value| cfg.dt = value;
    if (shared.frames) |frames| {
        cfg.output_interval = common.framesToInterval(cfg.steps, frames);
    }
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
