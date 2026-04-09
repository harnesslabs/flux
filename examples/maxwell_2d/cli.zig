//! `flux-examples maxwell-2d` subcommand: 2D electromagnetic simulation on
//! simplicial meshes. Supports a dipole-driven cavity demo and an exact TE₁₀
//! standing-wave demo.

const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const maxwell = @import("example_physics");

const Mesh2D = flux.topology.Mesh(2, 2);
const MaxwellState = maxwell.State(Mesh2D);
const PointDipole = maxwell.PointDipole(Mesh2D);

const flux_io = flux.io;

// ═══════════════════════════════════════════════════════════════════════════
// Configuration
// ═══════════════════════════════════════════════════════════════════════════

const Demo = enum { dipole, cavity };

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

// ═══════════════════════════════════════════════════════════════════════════
// Entry point
// ═══════════════════════════════════════════════════════════════════════════

pub fn run(allocator: std.mem.Allocator, args: []const [:0]const u8) !void {
    var config = Config{};
    var co = common.Common{};
    var parser = common.Parser.init(args);
    _ = parser.next();

    while (parser.next()) |arg| {
        if (eql(arg, "--demo")) {
            const val = try parser.requireValue("--demo");
            if (eql(val, "dipole")) {
                config.demo = .dipole;
            } else if (eql(val, "cavity")) {
                config.demo = .cavity;
            } else {
                std.debug.print("error: unknown demo '{s}'. available: dipole, cavity\n", .{val});
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
        if (try parser.tryCommon(arg, &co)) continue;
        std.debug.print("error: unknown flag '{s}'\n\n", .{arg});
        printUsage();
        std.process.exit(2);
    }

    if (co.help) {
        printUsage();
        return;
    }

    applyCommon(&config, co);

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const stderr = &stderr_writer.interface;

    switch (config.demo) {
        .dipole => try runDipole(allocator, config, stderr),
        .cavity => try runCavity(allocator, config, stderr),
    }
}

fn applyCommon(cfg: *Config, co: common.Common) void {
    common.applySharedFields(cfg, co);
    if (co.dt) |dt_value| {
        // 2D Maxwell derives dt from courant * h. Solve for the courant
        // number that produces the requested dt so dt() stays the single
        // source of truth.
        const h = cfg.domain / @as(f64, @floatFromInt(cfg.grid));
        cfg.courant = dt_value / h;
    }
}

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
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

    var series = try common.Series.init(
        allocator,
        config.output_dir,
        base_name,
        common.Plan.fromFrames(config.steps, config.frames, .{}),
    );
    defer series.deinit();

    var progress = common.Progress(@TypeOf(writer)).init(writer, config.steps);

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
        if (series.dueAt(@intCast(step_idx + 1), config.steps)) {
            try series.capture(t + time_step, Maxwell2DRenderer{ .state = state });
        }

        // Safe: step_idx < config.steps (u32), so step_idx + 1 ≤ u32 max.
        progress.update(@intCast(step_idx + 1));
    }

    progress.finish();
    const elapsed_s = progress.elapsed();

    try series.finalize();

    const energy_final = try maxwell.electromagnetic_energy(allocator, state);

    return .{
        .elapsed_s = elapsed_s,
        .energy_final = energy_final,
        .snapshot_count = series.count,
    };
}

const Maxwell2DRenderer = struct {
    state: *const MaxwellState,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try flux_io.write_fields(
            allocator,
            writer,
            self.state.mesh.*,
            self.state.E.values,
            self.state.B.values,
        );
    }
};

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

    const result = try simulate(allocator, &state, dipole, config, "dipole", writer);
    const t_final = @as(f64, @floatFromInt(config.steps)) * config.dt();
    const steps_per_sec = @as(f64, @floatFromInt(config.steps)) / result.elapsed_s;

    try writer.writeAll("\n  ── Results ─────────────────────────────────\n\n");
    try writer.print("  time     {d:.4}\n", .{t_final});
    try writer.print("  energy   {d:.6}\n", .{result.energy_final});
    var dur_buf: [16]u8 = undefined;
    try writer.print("  elapsed  {s} ({d:.0} steps/s)\n", .{ common.formatDuration(&dur_buf, result.elapsed_s), steps_per_sec });
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

    const result = try simulate(allocator, &state, null, config, "cavity", writer);
    const t_final = @as(f64, @floatFromInt(config.steps)) * config.dt();
    const periods = t_final / field_period;
    const drift_pct = @abs(result.energy_final - energy_initial) / energy_initial * 100.0;
    const steps_per_sec = @as(f64, @floatFromInt(config.steps)) / result.elapsed_s;

    try writer.writeAll("\n  ── Results ─────────────────────────────────\n\n");
    try writer.print("  time     {d:.4} ({d:.2} periods)\n", .{ t_final, periods });
    try writer.print("  energy   {d:.6} (drift {d:.4}%)\n", .{ result.energy_final, drift_pct });
    var dur_buf: [16]u8 = undefined;
    try writer.print("  elapsed  {s} ({d:.0} steps/s)\n", .{ common.formatDuration(&dur_buf, result.elapsed_s), steps_per_sec });
    if (result.snapshot_count > 0) {
        try writer.print("  output   {d} frames → {s}/\n", .{ result.snapshot_count, config.output_dir });
        try writer.print("\n  ▸ uv run tools/visualize.py {s}\n\n", .{config.output_dir});
    } else {
        try writer.writeAll("\n");
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Usage
// ═══════════════════════════════════════════════════════════════════════════

fn printUsage() void {
    std.debug.print(
        \\
        \\  flux-examples maxwell-2d — 2D electromagnetic simulation on simplicial meshes
        \\
        \\  usage: flux-examples maxwell-2d [--demo <name>] [options]
        \\
        \\  demos:
        \\    dipole    (default) point dipole radiating in a PEC cavity
        \\    cavity    TE₁₀ standing wave — exact mode, source-free
        \\
        \\  mesh & physics:
        \\    --grid N          cells per side (default: 32)
        \\    --domain L        domain side length (default: 1.0)
        \\    --courant C       Courant number, dt = C·h (default: 0.1)
        \\    --dt DT           explicit timestep override (computes courant = dt/h)
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
    , .{});
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
