//! flux CLI entry point.
//!
//! Usage:
//!   zig build run -- [--steps N] [--output DIR]          # dipole simulation (default)
//!   zig build run -- --demo cavity [--steps N] [--output DIR]  # cavity resonance
//!
//! When no --demo is specified, runs a point dipole simulation on a 32×32
//! PEC cavity — the simplest way to see the Maxwell solver in action.
//!
//! Demos:
//!   dipole   (default) Point dipole radiating in a PEC cavity.
//!   cavity   TE₁₀ standing wave in a PEC resonant cavity.

const std = @import("std");
const flux = @import("flux");

const Mesh2D = flux.Mesh(2);
const MaxwellState = flux.MaxwellState(Mesh2D);
const MaxwellRunner = flux.Runner(Mesh2D);
const PointDipole = flux.PointDipole(Mesh2D);

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    // Parse arguments.
    var demo_name: ?[]const u8 = null;
    var steps: u32 = 1000;
    var output_dir: []const u8 = "output";

    var i: usize = 1; // skip argv[0]
    while (i < args.len) : (i += 1) {
        if (std.mem.eql(u8, args[i], "--demo") and i + 1 < args.len) {
            i += 1;
            demo_name = args[i];
        } else if (std.mem.eql(u8, args[i], "--steps") and i + 1 < args.len) {
            i += 1;
            steps = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("error: invalid --steps value: {s}\n", .{args[i]});
                return error.InvalidArgument;
            };
        } else if (std.mem.eql(u8, args[i], "--output") and i + 1 < args.len) {
            i += 1;
            output_dir = args[i];
        } else if (std.mem.eql(u8, args[i], "--help") or std.mem.eql(u8, args[i], "-h")) {
            printUsage();
            return;
        } else {
            std.debug.print("error: unknown argument: {s}\n", .{args[i]});
            printUsage();
            return error.InvalidArgument;
        }
    }

    // Default to dipole simulation when no --demo is specified.
    const selected_demo = demo_name orelse "dipole";

    if (std.mem.eql(u8, selected_demo, "dipole")) {
        try runDipoleDemo(allocator, steps, output_dir);
    } else if (std.mem.eql(u8, selected_demo, "cavity")) {
        try runCavityDemo(allocator, steps, output_dir);
    } else {
        std.debug.print("error: unknown demo: {s}\n", .{selected_demo});
        std.debug.print("available demos: dipole, cavity\n", .{});
        return error.InvalidArgument;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Dipole demo (default)
// ═══════════════════════════════════════════════════════════════════════════

/// Point dipole radiating inside a [0,1]² PEC cavity.
///
/// A sinusoidal current source at the cavity center drives electromagnetic
/// waves that reflect off the PEC walls. This is the default simulation —
/// the simplest way to see the Maxwell solver producing non-trivial fields.
fn runDipoleDemo(allocator: std.mem.Allocator, steps: u32, output_dir: []const u8) !void {
    const domain_length: f64 = 1.0;
    const grid_n: u32 = 32;
    const grid_spacing = domain_length / @as(f64, @floatFromInt(grid_n));

    // CFL-stable timestep (Courant number 0.1).
    const dt = 0.1 * grid_spacing;

    // Source frequency: TE₁₀ resonance ω = π/L, f = 1/(2L).
    const source_freq = 1.0 / (2.0 * domain_length);
    const source_amplitude: f64 = 1.0;

    std.debug.print("═══ Dipole Simulation ═══\n\n", .{});
    std.debug.print("  domain:     [0, {d:.1}] × [0, {d:.1}]\n", .{ domain_length, domain_length });
    std.debug.print("  grid:       {d}×{d} ({d} cells, h = {d:.6})\n", .{ grid_n, grid_n, 2 * grid_n * grid_n, grid_spacing });
    std.debug.print("  dt:         {d:.6} (Courant = 0.1)\n", .{dt});
    std.debug.print("  steps:      {d}\n", .{steps});
    std.debug.print("  source:     dipole at center, f = {d:.4} Hz, A = {d:.1}\n", .{ source_freq, source_amplitude });
    std.debug.print("  output:     {s}/\n\n", .{output_dir});

    // Build mesh.
    var mesh = try Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);

    std.debug.print("  mesh:       {d} vertices, {d} edges, {d} faces\n", .{
        mesh.num_vertices(),
        mesh.num_edges(),
        mesh.num_faces(),
    });

    // Initialize zero fields.
    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Place dipole at cavity center.
    const center = [2]f64{ domain_length / 2.0, domain_length / 2.0 };
    const dipole = PointDipole.init(&mesh, source_freq, source_amplitude, center);

    std.debug.print("  dipole:     edge {d} (length = {d:.6})\n\n", .{
        dipole.edge_index,
        dipole.edge_length,
    });

    // Ensure output directory exists.
    std.fs.cwd().makeDir(output_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    // Write VTK snapshots every 10 steps.
    const output_interval: u32 = @max(1, steps / 100);

    std.debug.print("  simulating {d} steps...\n", .{steps});

    // Run simulation.
    try MaxwellRunner.run(allocator, &state, dipole, .{
        .steps = steps,
        .dt = dt,
        .output_interval = output_interval,
        .output_path = output_dir,
        .output_base_name = "dipole",
    });

    // Final diagnostics.
    const energy_final = try flux.electromagnetic_energy(allocator, &state);
    const t_final = @as(f64, @floatFromInt(steps)) * dt;

    std.debug.print("\n  ── results ──\n", .{});
    std.debug.print("  t_final:    {d:.4}\n", .{t_final});
    std.debug.print("  energy(T):  {d:.6}\n", .{energy_final});
    std.debug.print("  snapshots:  {s}/dipole_*.vtu\n", .{output_dir});
    std.debug.print("  animation:  {s}/dipole.pvd\n\n", .{output_dir});
    std.debug.print("  Open {s}/dipole.pvd in ParaView to view the radiation pattern.\n", .{output_dir});
}

// ═══════════════════════════════════════════════════════════════════════════
// Cavity demo
// ═══════════════════════════════════════════════════════════════════════════

/// TE₁₀ standing wave in a [0,1]² PEC resonant cavity.
///
/// Analytical solution (c = 1):
///   E_y(x,y,t) = sin(πx) sin(πt),  E_x = 0
///   B_z(x,y,t) = cos(πx) cos(πt)
///   ω = π,  f = 1/2,  T = 2
///
/// The demo initializes the exact mode and runs source-free with PEC
/// boundary conditions. VTK snapshots show the standing wave pattern.
fn runCavityDemo(allocator: std.mem.Allocator, steps: u32, output_dir: []const u8) !void {
    const domain_length: f64 = 1.0;
    const grid_n: u32 = 32;
    const grid_spacing = domain_length / @as(f64, @floatFromInt(grid_n));

    // CFL-stable timestep (Courant number 0.1).
    const dt = 0.1 * grid_spacing;
    const field_period = 2.0 * domain_length; // T = 2L/c with c = 1
    const analytical_freq = 1.0 / field_period;

    std.debug.print("═══ TE₁₀ Cavity Resonance Demo ═══\n\n", .{});
    std.debug.print("  domain:     [0, {d:.1}] × [0, {d:.1}]\n", .{ domain_length, domain_length });
    std.debug.print("  grid:       {d}×{d} ({d} cells, h = {d:.6})\n", .{ grid_n, grid_n, 2 * grid_n * grid_n, grid_spacing });
    std.debug.print("  dt:         {d:.6} (Courant = 0.1)\n", .{dt});
    std.debug.print("  steps:      {d}\n", .{steps});
    std.debug.print("  period:     T = {d:.4} (f = {d:.4}, ω = π/{d:.1})\n", .{ field_period, analytical_freq, domain_length });
    std.debug.print("  output:     {s}/\n\n", .{output_dir});

    // Build mesh.
    var mesh = try Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);

    std.debug.print("  mesh:       {d} vertices, {d} edges, {d} faces\n", .{
        mesh.num_vertices(),
        mesh.num_edges(),
        mesh.num_faces(),
    });

    // Initialize Maxwell state with TE₁₀ mode.
    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // B at t = -dt/2 (leapfrog stagger), E at t = 0 (zero = sin(0)).
    flux.project_te10_b(&mesh, state.B.values, -dt / 2.0, domain_length);

    // Compute initial energy.
    const energy_initial = try flux.electromagnetic_energy(allocator, &state);
    std.debug.print("  energy(0):  {d:.6}\n\n", .{energy_initial});

    // Ensure output directory exists.
    std.fs.cwd().makeDir(output_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    // Write VTK snapshots every ~T/20 (20 frames per period).
    const frames_per_period: u32 = 20;
    const steps_per_period: u32 = @intFromFloat(@round(field_period / dt));
    const output_interval = @max(1, steps_per_period / frames_per_period);

    std.debug.print("  simulating {d} steps...\n", .{steps});

    // Run simulation with VTK output.
    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = steps,
        .dt = dt,
        .output_interval = output_interval,
        .output_path = output_dir,
        .output_base_name = "cavity",
    });

    // Final diagnostics.
    const energy_final = try flux.electromagnetic_energy(allocator, &state);
    const energy_ratio = energy_final / energy_initial;
    const t_final = @as(f64, @floatFromInt(steps)) * dt;
    const periods_simulated = t_final / field_period;

    std.debug.print("\n  ── results ──\n", .{});
    std.debug.print("  t_final:    {d:.4} ({d:.2} periods)\n", .{ t_final, periods_simulated });
    std.debug.print("  energy(T):  {d:.6} (ratio = {d:.6})\n", .{ energy_final, energy_ratio });
    std.debug.print("  snapshots:  {s}/cavity_*.vtu\n", .{output_dir});
    std.debug.print("  animation:  {s}/cavity.pvd\n\n", .{output_dir});
    std.debug.print("  Open {s}/cavity.pvd in ParaView to view the standing wave.\n", .{output_dir});
}

fn printUsage() void {
    std.debug.print(
        \\usage: flux [--demo <name>] [options]
        \\
        \\demos:
        \\  dipole    (default) point dipole radiating in a PEC cavity
        \\  cavity    TE₁₀ standing wave in a PEC resonant cavity
        \\
        \\options:
        \\  --steps N     number of timesteps (default: 1000)
        \\  --output DIR  output directory for VTK files (default: output)
        \\  --help        show this message
        \\
    , .{});
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
