//! maxwell_2d — 2D electromagnetic simulation on simplicial meshes.
//!
//! Runs electromagnetic simulations on a 2D triangulated PEC cavity and writes
//! VTK snapshots for visualization. Every physical and numerical parameter is
//! configurable from the command line.
//!
//! Usage:
//!   zig build example-maxwell2d -- [options]                     # dipole (default)
//!   zig build example-maxwell2d -- --demo cavity [options]       # standing wave
//!   zig build example-maxwell2d -- --help                        # full flag reference

const std = @import("std");
const flux = @import("flux");
const maxwell = @import("maxwell.zig");

const Mesh2D = flux.Mesh(2);
const MaxwellState = maxwell.State(Mesh2D);
const PointDipole = maxwell.PointDipole(Mesh2D);

const flux_io = flux.io;
const std_fs = std.fs;

// ═══════════════════════════════════════════════════════════════════════════
// Configuration
// ═══════════════════════════════════════════════════════════════════════════

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

    fn spacing(self: Config) f64 {
        return self.domain / @as(f64, @floatFromInt(self.grid));
    }

    fn dt(self: Config) f64 {
        return self.courant * self.spacing();
    }

    fn outputInterval(self: Config) u32 {
        if (self.frames == 0) return 0;
        return @max(1, self.steps / self.frames);
    }

    fn sourceFrequency(self: Config) f64 {
        if (self.frequency != 0.0) return self.frequency;
        return 1.0 / (2.0 * self.domain); // TE₁₀ resonance
    }
};

const Demo = enum { dipole, cavity };

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
    const stderr = (std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter();

    switch (config.demo) {
        .dipole => try runDipole(allocator, config, stderr),
        .cavity => try runCavity(allocator, config, stderr),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Progress display
// ═══════════════════════════════════════════════════════════════════════════

fn Progress(comptime Writer: type) type {
    return struct {
        const Self = @This();

        writer: Writer,
        total: u32,
        timer: std.time.Timer,
        last_draw_ns: u64 = 0,
        bar_width: u32 = 40,

        fn init(writer: Writer, total: u32) Self {
            return .{
                .writer = writer,
                .total = total,
                .timer = std.time.Timer.start() catch
                    @panic("OS timer unavailable — cannot run simulation"),
            };
        }

        fn update(self: *Self, step: u32) void {
            const elapsed_ns = self.timer.read();

            // Avoid flooding stderr on fast machines — cap at ~20 fps.
            if (elapsed_ns - self.last_draw_ns < 50_000_000 and step < self.total) return;
            self.last_draw_ns = elapsed_ns;

            const elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0;
            const frac = @as(f64, @floatFromInt(step)) / @as(f64, @floatFromInt(self.total));
            const pct = frac * 100.0;

            // Guard against division by zero in the first few milliseconds.
            const steps_per_sec = if (elapsed_s > 0.01) @as(f64, @floatFromInt(step)) / elapsed_s else 0.0;
            const remaining = @as(f64, @floatFromInt(self.total - step));
            const eta_s = if (steps_per_sec > 0.01) remaining / steps_per_sec else 0.0;

            const filled: u32 = @intFromFloat(frac * @as(f64, @floatFromInt(self.bar_width)));

            // Stack buffer avoids allocation on every redraw.
            var bar: [64]u8 = undefined;
            for (0..self.bar_width) |j| {
                bar[j] = if (j < filled) '#' else '-';
            }
            const bar_str = bar[0..self.bar_width];

            var elapsed_buf: [16]u8 = undefined;
            var eta_buf: [16]u8 = undefined;
            self.writer.print("\r  {s}  {d:>5.1}%  {d}/{d}  {s}  ETA {s}  {d:.0} steps/s    ", .{
                bar_str,
                pct,
                step,
                self.total,
                formatDuration(&elapsed_buf, elapsed_s),
                formatDuration(&eta_buf, eta_s),
                steps_per_sec,
            }) catch return;
        }

        fn finish(self: *Self) void {
            // Clear the progress line.
            self.writer.writeAll("\r") catch return;
            for (0..120) |_| self.writer.writeByte(' ') catch return;
            self.writer.writeAll("\r") catch return;
        }

        fn elapsed(self: *Self) f64 {
            return @as(f64, @floatFromInt(self.timer.read())) / 1_000_000_000.0;
        }
    };
}

fn formatDuration(buf: *[16]u8, seconds: f64) []const u8 {
    if (seconds < 60.0) {
        return std.fmt.bufPrint(buf, "{d:.1}s", .{seconds}) catch "??";
    }
    const mins: u32 = @intFromFloat(seconds / 60.0);
    const secs: u32 = @intFromFloat(@mod(seconds, 60.0));
    return std.fmt.bufPrint(buf, "{d}m{d:0>2}s", .{ mins, secs }) catch "??";
}

// ═══════════════════════════════════════════════════════════════════════════
// VTK snapshot helpers
// ═══════════════════════════════════════════════════════════════════════════

fn writeSnapshot(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    filename: []const u8,
    state: *const MaxwellState,
) !void {
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try flux_io.write_fields(
        allocator,
        output.writer(allocator),
        Mesh2D.dimension,
        state.mesh.*,
        state.E.values,
        state.B.values,
    );

    var dir = try std_fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(filename, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn writePvd(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    base_name: []const u8,
    entries: []const flux_io.PvdEntry,
) !void {
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try flux_io.write_pvd(output.writer(allocator), entries);

    var pvd_buf: [flux_io.max_snapshot_filename_length]u8 = undefined;
    const pvd_name = std.fmt.bufPrint(&pvd_buf, "{s}.pvd", .{base_name}) catch
        return error.FilenameTooLong;

    var dir = try std_fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(pvd_name, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn ensureDir(path: []const u8) !void {
    std.fs.cwd().makeDir(path) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Simulation loop (shared by both demos)
// ═══════════════════════════════════════════════════════════════════════════

const SimResult = struct {
    elapsed_s: f64,
    energy_final: f64,
    snapshot_count: u32,
};

fn simulate(
    allocator: std.mem.Allocator,
    state: *MaxwellState,
    source: ?PointDipole,
    config: Config,
    base_name: []const u8,
    writer: anytype,
) !SimResult {
    const time_step = config.dt();
    const interval = config.outputInterval();
    const has_output = interval > 0;

    // PVD entry storage. Saturating add guards the pathological case where
    // steps == u32_max and interval == 1.
    const max_snapshots: u32 = if (has_output) (config.steps / interval) +| 1 else 0;
    var pvd_entries: []flux_io.PvdEntry = &.{};
    var filename_bufs: [][flux_io.max_snapshot_filename_length]u8 = &.{};
    var snapshot_count: u32 = 0;

    if (has_output) {
        pvd_entries = try allocator.alloc(flux_io.PvdEntry, max_snapshots);
        filename_bufs = try allocator.alloc([flux_io.max_snapshot_filename_length]u8, max_snapshots);
    }
    defer if (has_output) {
        allocator.free(pvd_entries);
        allocator.free(filename_bufs);
    };

    var progress = Progress(@TypeOf(writer)).init(writer, config.steps);

    for (0..config.steps) |step_idx| {
        const t = @as(f64, @floatFromInt(step_idx)) * time_step;

        // Source.
        if (source) |s| {
            s.apply(&state.J, t);
        }

        // Advance fields.
        try maxwell.leapfrog_step(allocator, state, time_step);
        maxwell.apply_pec_boundary(state);

        // VTK snapshot.
        if (has_output and (step_idx + 1) % interval == 0) {
            const filename = flux_io.snapshot_filename(
                &filename_bufs[snapshot_count],
                base_name,
                snapshot_count,
            );
            try writeSnapshot(allocator, config.output_dir, filename, state);
            pvd_entries[snapshot_count] = .{
                .timestep = t + time_step,
                .filename = filename,
            };
            snapshot_count += 1;
        }

        // Safe: step_idx < config.steps (u32), so step_idx + 1 ≤ u32 max.
        progress.update(@intCast(step_idx + 1));
    }

    progress.finish();
    const elapsed_s = progress.elapsed();

    // Write PVD collection.
    if (has_output and snapshot_count > 0) {
        try writePvd(allocator, config.output_dir, base_name, pvd_entries[0..snapshot_count]);
    }

    const energy_final = try maxwell.electromagnetic_energy(allocator, state);

    return .{
        .elapsed_s = elapsed_s,
        .energy_final = energy_final,
        .snapshot_count = snapshot_count,
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Dipole demo
// ═══════════════════════════════════════════════════════════════════════════

fn runDipole(allocator: std.mem.Allocator, config: Config, writer: anytype) !void {
    const h = config.spacing();
    const time_step = config.dt();
    const freq = config.sourceFrequency();

    try writer.writeAll("\n  ── Dipole Simulation ───────────────────────\n\n");
    try writer.print("  domain    [0, {d:.2}]²\n", .{config.domain});
    try writer.print("  grid      {d}×{d} ({d} triangles)\n", .{ config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid });
    try writer.print("  spacing   h = {d:.6}\n", .{h});
    try writer.print("  timestep  dt = {d:.6} (Courant {d:.2})\n", .{ time_step, config.courant });
    try writer.print("  source    f = {d:.4} Hz, A = {d:.2}\n\n", .{ freq, config.amplitude });

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    try writer.print("  mesh     {d} vertices  {d} edges  {d} faces\n", .{
        mesh.num_vertices(), mesh.num_edges(), mesh.num_faces(),
    });

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    const center = [2]f64{ config.domain / 2.0, config.domain / 2.0 };
    const dipole = PointDipole.init(&mesh, freq, config.amplitude, center);
    try writer.print("  source   edge {d} (length {d:.6})\n\n", .{ dipole.edge_index, dipole.edge_length });

    try ensureDir(config.output_dir);

    const result = try simulate(allocator, &state, dipole, config, "dipole", writer);
    const t_final = @as(f64, @floatFromInt(config.steps)) * config.dt();
    const steps_per_sec = @as(f64, @floatFromInt(config.steps)) / result.elapsed_s;

    try writer.writeAll("\n  ── Results ─────────────────────────────────\n\n");
    try writer.print("  time     {d:.4}\n", .{t_final});
    try writer.print("  energy   {d:.6}\n", .{result.energy_final});
    var dur_buf: [16]u8 = undefined;
    try writer.print("  elapsed  {s} ({d:.0} steps/s)\n", .{ formatDuration(&dur_buf, result.elapsed_s), steps_per_sec });
    if (result.snapshot_count > 0) {
        try writer.print("  output   {d} frames → {s}/\n", .{ result.snapshot_count, config.output_dir });
        try writer.print("\n  ▸ uv run tools/visualize.py {s}\n\n", .{config.output_dir});
    } else {
        try writer.writeAll("\n");
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Cavity demo
// ═══════════════════════════════════════════════════════════════════════════

fn runCavity(allocator: std.mem.Allocator, config: Config, writer: anytype) !void {
    const h = config.spacing();
    const time_step = config.dt();
    const field_period = 2.0 * config.domain;
    const omega = std.math.pi / config.domain;

    try writer.writeAll("\n  ── TE₁₀ Cavity Resonance ──────────────────\n\n");
    try writer.print("  domain    [0, {d:.2}]²\n", .{config.domain});
    try writer.print("  grid      {d}×{d} ({d} triangles)\n", .{ config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid });
    try writer.print("  spacing   h = {d:.6}\n", .{h});
    try writer.print("  timestep  dt = {d:.6} (Courant {d:.2})\n", .{ time_step, config.courant });
    try writer.print("  period    T = {d:.4} (ω = {d:.4})\n\n", .{ field_period, omega });

    var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    try writer.print("  mesh     {d} vertices  {d} edges  {d} faces\n", .{
        mesh.num_vertices(), mesh.num_edges(), mesh.num_faces(),
    });

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // B at t = -dt/2 (leapfrog stagger), E at t = 0 is zero.
    maxwell.project_te10_b(&mesh, state.B.values, -time_step / 2.0, config.domain);

    const energy_initial = try maxwell.electromagnetic_energy(allocator, &state);
    try writer.print("  energy₀  {d:.6}\n\n", .{energy_initial});

    try ensureDir(config.output_dir);

    const result = try simulate(allocator, &state, null, config, "cavity", writer);
    const t_final = @as(f64, @floatFromInt(config.steps)) * config.dt();
    const periods = t_final / field_period;
    const drift_pct = @abs(result.energy_final - energy_initial) / energy_initial * 100.0;
    const steps_per_sec = @as(f64, @floatFromInt(config.steps)) / result.elapsed_s;

    try writer.writeAll("\n  ── Results ─────────────────────────────────\n\n");
    try writer.print("  time     {d:.4} ({d:.2} periods)\n", .{ t_final, periods });
    try writer.print("  energy   {d:.6} (drift {d:.4}%)\n", .{ result.energy_final, drift_pct });
    var dur_buf: [16]u8 = undefined;
    try writer.print("  elapsed  {s} ({d:.0} steps/s)\n", .{ formatDuration(&dur_buf, result.elapsed_s), steps_per_sec });
    if (result.snapshot_count > 0) {
        try writer.print("  output   {d} frames → {s}/\n", .{ result.snapshot_count, config.output_dir });
        try writer.print("\n  ▸ uv run tools/visualize.py {s}\n\n", .{config.output_dir});
    } else {
        try writer.writeAll("\n");
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Argument parsing
// ═══════════════════════════════════════════════════════════════════════════

const ParseError = error{InvalidArgument};

fn parseArgs(args: []const [:0]const u8) ParseError!Config {
    var config = Config{};
    var i: usize = 1;
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
            std.debug.print("error: unknown argument '{s}'\n\n", .{arg});
            printUsage();
            return ParseError.InvalidArgument;
        }
    }
    return config;
}

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
// Usage
// ═══════════════════════════════════════════════════════════════════════════

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux — 2D electromagnetic simulation on simplicial meshes
        \\
        \\  usage: flux [--demo <name>] [options]
        \\
        \\  demos:
        \\    dipole    (default) point dipole radiating in a PEC cavity
        \\    cavity    TE₁₀ standing wave — exact mode, source-free
        \\
        \\  mesh & physics:
        \\    --grid N          cells per side (default: 32)
        \\    --domain L        domain side length (default: 1.0)
        \\    --courant C       Courant number, dt = C·h (default: 0.1)
        \\
        \\  dipole source (ignored for cavity):
        \\    --frequency F     source frequency in Hz (default: TE₁₀ resonance)
        \\    --amplitude A     source amplitude (default: 1.0)
        \\
        \\  time stepping:
        \\    --steps N         number of timesteps (default: 1000)
        \\
        \\  output:
        \\    --output DIR      VTK output directory (default: output)
        \\    --frames N        number of snapshots (default: 100)
        \\
        \\  examples:
        \\    zig build run                                    # defaults
        \\    zig build run -- --demo cavity --steps 2000      # standing wave
        \\    zig build run -- --grid 64 --steps 4000          # fine mesh
        \\    zig build run -- --frequency 1.5 --amplitude 2   # off-resonance
        \\    zig build run -- --grid 16 --courant 0.2         # coarse & fast
        \\
        \\  visualization:
        \\    uv run tools/visualize.py output --field B_flux --output anim.gif
        \\
    , .{});
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
    // Pull in maxwell module tests (physics integration tests).
    _ = @import("maxwell.zig");
}
