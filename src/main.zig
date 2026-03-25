//! flux CLI entry point.
//!
//! Runs electromagnetic simulations on a 2D triangulated PEC cavity and writes
//! VTK snapshots for visualization. Every physical and numerical parameter is
//! configurable from the command line.
//!
//! Usage:
//!   zig build run -- [options]                     # dipole (default)
//!   zig build run -- --demo cavity [options]       # standing wave
//!   zig build run -- --help                        # full flag reference

const std = @import("std");
const flux = @import("flux");

const Mesh2D = flux.Mesh(2);
const MaxwellState = flux.MaxwellState(Mesh2D);
const MaxwellRunner = flux.Runner(Mesh2D);
const PointDipole = flux.PointDipole(Mesh2D);

// ═══════════════════════════════════════════════════════════════════════════
// Configuration
// ═══════════════════════════════════════════════════════════════════════════

/// All simulation parameters, parsed from CLI args with sensible defaults.
const Config = struct {
    demo: Demo = .dipole,
    steps: u32 = 1000,
    grid: u32 = 32,
    domain: f64 = 1.0,
    courant: f64 = 0.1,
    frequency: f64 = 0.0, // 0 = auto (TE₁₀ resonance)
    amplitude: f64 = 1.0,
    output_dir: []const u8 = "output",
    frames: u32 = 100,
};

const Demo = enum { dipole, cavity };

const ParseError = error{InvalidArgument};

fn parseArgs(args: []const [:0]const u8) ParseError!Config {
    var config = Config{};
    var i: usize = 1; // skip argv[0]
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (eql(arg, "--help") or eql(arg, "-h")) {
            printUsage();
            std.process.exit(0);
        } else if (eql(arg, "--demo")) {
            const val = nextArg(args, &i) orelse return flagError("--demo");
            if (eql(val, "dipole")) {
                config.demo = .dipole;
            } else if (eql(val, "cavity")) {
                config.demo = .cavity;
            } else {
                std.debug.print("error: unknown demo '{s}'. available: dipole, cavity\n", .{val});
                return ParseError.InvalidArgument;
            }
        } else if (eql(arg, "--steps")) {
            config.steps = parseU32(args, &i, "--steps") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--grid")) {
            config.grid = parseU32(args, &i, "--grid") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--domain")) {
            config.domain = parseF64(args, &i, "--domain") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--courant")) {
            config.courant = parseF64(args, &i, "--courant") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--frequency")) {
            config.frequency = parseF64(args, &i, "--frequency") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--amplitude")) {
            config.amplitude = parseF64(args, &i, "--amplitude") orelse return ParseError.InvalidArgument;
        } else if (eql(arg, "--output")) {
            config.output_dir = nextArg(args, &i) orelse return flagError("--output");
        } else if (eql(arg, "--frames")) {
            config.frames = parseU32(args, &i, "--frames") orelse return ParseError.InvalidArgument;
        } else {
            std.debug.print("error: unknown argument '{s}'\n", .{arg});
            printUsage();
            return ParseError.InvalidArgument;
        }
    }
    return config;
}

// ── Argument parsing helpers ─────────────────────────────────────────────

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
    const val = nextArg(args, i) orelse {
        std.debug.print("error: {s} requires a value\n", .{flag});
        return null;
    };
    return std.fmt.parseInt(u32, val, 10) catch {
        std.debug.print("error: invalid {s} value: {s}\n", .{ flag, val });
        return null;
    };
}

fn parseF64(args: []const [:0]const u8, i: *usize, flag: []const u8) ?f64 {
    const val = nextArg(args, i) orelse {
        std.debug.print("error: {s} requires a value\n", .{flag});
        return null;
    };
    return std.fmt.parseFloat(f64, val) catch {
        std.debug.print("error: invalid {s} value: {s}\n", .{ flag, val });
        return null;
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Entry point
// ═══════════════════════════════════════════════════════════════════════════

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    const config = parseArgs(args) catch return;

    switch (config.demo) {
        .dipole => try runDipoleDemo(allocator, config),
        .cavity => try runCavityDemo(allocator, config),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Shared helpers
// ═══════════════════════════════════════════════════════════════════════════

fn ensureDir(path: []const u8) !void {
    std.fs.cwd().makeDir(path) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };
}

fn computeDt(config: Config) f64 {
    const h = config.domain / @as(f64, @floatFromInt(config.grid));
    return config.courant * h;
}

fn outputInterval(config: Config) u32 {
    if (config.frames == 0) return 0;
    return @max(1, config.steps / config.frames);
}

fn printMeshInfo(mesh: *const Mesh2D) void {
    std.debug.print("  mesh:       {d} vertices, {d} edges, {d} faces\n", .{
        mesh.num_vertices(),
        mesh.num_edges(),
        mesh.num_faces(),
    });
}

// ═══════════════════════════════════════════════════════════════════════════
// Dipole demo
// ═══════════════════════════════════════════════════════════════════════════

fn runDipoleDemo(allocator: std.mem.Allocator, config: Config) !void {
    const h = config.domain / @as(f64, @floatFromInt(config.grid));
    const dt = computeDt(config);

    // Default frequency: TE₁₀ resonance f = 1/(2L).
    const freq = if (config.frequency == 0.0) 1.0 / (2.0 * config.domain) else config.frequency;

    std.debug.print("═══ Dipole Simulation ═══\n\n", .{});
    std.debug.print("  domain:     [0, {d:.2}] × [0, {d:.2}]\n", .{ config.domain, config.domain });
    std.debug.print("  grid:       {d}×{d} ({d} cells, h = {d:.6})\n", .{ config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid, h });
    std.debug.print("  dt:         {d:.6} (Courant = {d:.2})\n", .{ dt, config.courant });
    std.debug.print("  steps:      {d}\n", .{config.steps});
    std.debug.print("  source:     f = {d:.4} Hz, A = {d:.2}\n", .{ freq, config.amplitude });
    std.debug.print("  output:     {s}/\n\n", .{config.output_dir});

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);
    printMeshInfo(&mesh);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    const center = [2]f64{ config.domain / 2.0, config.domain / 2.0 };
    const dipole = PointDipole.init(&mesh, freq, config.amplitude, center);

    std.debug.print("  dipole:     edge {d} (length = {d:.6})\n\n", .{
        dipole.edge_index,
        dipole.edge_length,
    });

    try ensureDir(config.output_dir);

    const interval = outputInterval(config);
    std.debug.print("  simulating {d} steps", .{config.steps});
    if (interval > 0) {
        std.debug.print(" (snapshot every {d})", .{interval});
    }
    std.debug.print("...\n", .{});

    try MaxwellRunner.run(allocator, &state, dipole, .{
        .steps = config.steps,
        .dt = dt,
        .output_interval = interval,
        .output_path = config.output_dir,
        .output_base_name = "dipole",
    });

    const energy_final = try flux.electromagnetic_energy(allocator, &state);
    const t_final = @as(f64, @floatFromInt(config.steps)) * dt;

    std.debug.print("\n  ── results ──\n", .{});
    std.debug.print("  t_final:    {d:.4}\n", .{t_final});
    std.debug.print("  energy(T):  {d:.6}\n", .{energy_final});
    std.debug.print("  snapshots:  {s}/dipole_*.vtu\n", .{config.output_dir});
    std.debug.print("  animation:  {s}/dipole.pvd\n\n", .{config.output_dir});
    std.debug.print("  Visualize:  uv run tools/visualize.py {s}\n", .{config.output_dir});
}

// ═══════════════════════════════════════════════════════════════════════════
// Cavity demo
// ═══════════════════════════════════════════════════════════════════════════

fn runCavityDemo(allocator: std.mem.Allocator, config: Config) !void {
    const h = config.domain / @as(f64, @floatFromInt(config.grid));
    const dt = computeDt(config);
    const field_period = 2.0 * config.domain; // T = 2L/c with c = 1
    const analytical_freq = 1.0 / field_period;

    std.debug.print("═══ TE₁₀ Cavity Resonance Demo ═══\n\n", .{});
    std.debug.print("  domain:     [0, {d:.2}] × [0, {d:.2}]\n", .{ config.domain, config.domain });
    std.debug.print("  grid:       {d}×{d} ({d} cells, h = {d:.6})\n", .{ config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid, h });
    std.debug.print("  dt:         {d:.6} (Courant = {d:.2})\n", .{ dt, config.courant });
    std.debug.print("  steps:      {d}\n", .{config.steps});
    std.debug.print("  period:     T = {d:.4} (f = {d:.4}, ω = π/{d:.2})\n", .{ field_period, analytical_freq, config.domain });
    std.debug.print("  output:     {s}/\n\n", .{config.output_dir});

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);
    printMeshInfo(&mesh);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // B at t = -dt/2 (leapfrog stagger), E at t = 0 (zero = sin(0)).
    flux.project_te10_b(&mesh, state.B.values, -dt / 2.0, config.domain);

    const energy_initial = try flux.electromagnetic_energy(allocator, &state);
    std.debug.print("  energy(0):  {d:.6}\n\n", .{energy_initial});

    try ensureDir(config.output_dir);

    const interval = outputInterval(config);
    std.debug.print("  simulating {d} steps", .{config.steps});
    if (interval > 0) {
        std.debug.print(" (snapshot every {d})", .{interval});
    }
    std.debug.print("...\n", .{});

    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = config.steps,
        .dt = dt,
        .output_interval = interval,
        .output_path = config.output_dir,
        .output_base_name = "cavity",
    });

    const energy_final = try flux.electromagnetic_energy(allocator, &state);
    const energy_ratio = energy_final / energy_initial;
    const t_final = @as(f64, @floatFromInt(config.steps)) * dt;
    const periods_simulated = t_final / field_period;

    std.debug.print("\n  ── results ──\n", .{});
    std.debug.print("  t_final:    {d:.4} ({d:.2} periods)\n", .{ t_final, periods_simulated });
    std.debug.print("  energy(T):  {d:.6} (ratio = {d:.6})\n", .{ energy_final, energy_ratio });
    std.debug.print("  snapshots:  {s}/cavity_*.vtu\n", .{config.output_dir});
    std.debug.print("  animation:  {s}/cavity.pvd\n\n", .{config.output_dir});
    std.debug.print("  Visualize:  uv run tools/visualize.py {s}\n", .{config.output_dir});
}

// ═══════════════════════════════════════════════════════════════════════════
// Usage
// ═══════════════════════════════════════════════════════════════════════════

fn printUsage() void {
    std.debug.print(
        \\usage: flux [--demo <name>] [options]
        \\
        \\Run a 2D electromagnetic simulation in a PEC cavity and write VTK
        \\output for visualization. All parameters have sensible defaults —
        \\just `zig build run` to get started.
        \\
        \\demos:
        \\  dipole    (default) point dipole radiating in a PEC cavity
        \\  cavity    TE₁₀ standing wave — exact analytical mode, source-free
        \\
        \\mesh & physics:
        \\  --grid N          grid cells per side (default: 32)
        \\  --domain L        square domain side length (default: 1.0)
        \\  --courant C       Courant number dt = C·h (default: 0.1)
        \\
        \\dipole source (ignored for cavity demo):
        \\  --frequency F     source frequency in Hz (default: TE₁₀ resonance)
        \\  --amplitude A     source amplitude (default: 1.0)
        \\
        \\time stepping:
        \\  --steps N         number of timesteps (default: 1000)
        \\
        \\output:
        \\  --output DIR      output directory for VTK files (default: output)
        \\  --frames N        number of VTK snapshots to write (default: 100)
        \\
        \\examples:
        \\  zig build run                                    # quick dipole, defaults
        \\  zig build run -- --demo cavity --steps 2000      # 3 periods of standing wave
        \\  zig build run -- --grid 64 --steps 4000          # finer mesh, longer run
        \\  zig build run -- --frequency 1.5 --amplitude 2   # off-resonance drive
        \\  zig build run -- --grid 16 --courant 0.2         # coarse & fast
        \\
        \\visualization:
        \\  uv run tools/visualize.py output --field B_flux --output anim.gif
        \\
    , .{});
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
