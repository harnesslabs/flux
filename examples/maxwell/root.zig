const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const runtime = @import("runtime.zig");
pub const reference = @import("reference.zig");

pub const BoundaryCondition = runtime.BoundaryCondition;

pub const DipoleConfig2D = struct {
    steps: u32 = 1000,
    grid: u32 = 32,
    domain: f64 = 1.0,
    courant: f64 = 0.1,
    frequency: f64 = 0.0,
    amplitude: f64 = 1.0,
    output_dir: []const u8 = "output",
    frames: u32 = 100,

    pub fn spacing(self: DipoleConfig2D) f64 {
        return self.domain / @as(f64, @floatFromInt(self.grid));
    }

    pub fn dt(self: DipoleConfig2D) f64 {
        return self.courant * self.spacing();
    }

    pub fn sourceFrequency(self: DipoleConfig2D) f64 {
        if (self.frequency != 0.0) return self.frequency;
        return 1.0 / (2.0 * self.domain);
    }
};

pub fn CavityConfig(comptime dim: u8) type {
    return switch (dim) {
        2 => struct {
            steps: u32 = 1000,
            grid: u32 = 32,
            domain: f64 = 1.0,
            courant: f64 = 0.1,
            output_dir: []const u8 = "output",
            frames: u32 = 100,

            pub fn spacing(self: @This()) f64 {
                return self.domain / @as(f64, @floatFromInt(self.grid));
            }

            pub fn dt(self: @This()) f64 {
                return self.courant * self.spacing();
            }
        },
        3 => struct {
            steps: u32 = 1000,
            nx: u32 = 2,
            ny: u32 = 2,
            nz: u32 = 2,
            width: f64 = 1.0,
            height: f64 = 1.0,
            depth: f64 = 1.0,
            dt: f64 = 0.01,
            output_dir: ?[]const u8 = null,
            output_interval: u32 = 0,

            pub fn gridSpacingMin(self: @This()) f64 {
                return @min(
                    self.width / @as(f64, @floatFromInt(self.nx)),
                    @min(
                        self.height / @as(f64, @floatFromInt(self.ny)),
                        self.depth / @as(f64, @floatFromInt(self.nz)),
                    ),
                );
            }
        },
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

const RunResult2D = struct {
    elapsed_s: f64,
    energy_final: f64,
    snapshot_count: u32,
};

const RunResult3D = struct {
    elapsed_s: f64,
    snapshot_count: u32,
};

pub fn Mesh(comptime dim: u8) type {
    return switch (dim) {
        2 => runtime.Mesh2D,
        3 => runtime.Mesh3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn System(comptime dim: u8) type {
    return runtime.Maxwell(dim, Mesh(dim));
}

pub fn RunResult(comptime dim: u8) type {
    return switch (dim) {
        2 => RunResult2D,
        3 => RunResult3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn runDipole(
    allocator: std.mem.Allocator,
    config: DipoleConfig2D,
    writer: *std.Io.Writer,
) !RunResult(2) {
    const SystemType = System(2);
    const dt = config.dt();
    var mesh = try runtime.Mesh2D.cartesian(
        allocator,
        .{ config.grid, config.grid },
        .{ config.domain, config.domain },
    );
    defer mesh.deinit(allocator);

    const frequency = config.sourceFrequency();
    try writer.writeAll("\n  ── Dipole Simulation ───────────────────────\n\n");
    try writer.print(
        "  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n  source    f = {d:.4} Hz, A = {d:.2}\n\n",
        .{
            config.domain,
            config.grid,
            config.grid,
            2 * @as(u64, config.grid) * config.grid,
            config.spacing(),
            dt,
            config.courant,
            frequency,
            config.amplitude,
        },
    );

    var system = try SystemType.dipole(allocator, &mesh, .{
        .center = .{ config.domain / 2.0, config.domain / 2.0 },
        .frequency = frequency,
        .amplitude = config.amplitude,
        .boundary = .pec,
    });
    defer system.deinit(allocator);

    const interval = @max(@as(u32, 1), common.framesToInterval(config.steps, config.frames));
    var evolution = try flux.evolution.Evolution(*SystemType, flux.integrators.Leapfrog).config()
        .dt(dt)
        .steps(config.steps)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(*SystemType)
                .field(.electric)
                .field(.magnetic)
                .measurement(.energy)
                .directory(config.output_dir)
                .baseName("dipole")
                .everySteps(interval),
        )
        .init(allocator, &system);
    defer evolution.deinit();

    const run_result = try evolution.run();
    return .{
        .elapsed_s = run_result.elapsed_s,
        .energy_final = try system.measurement(allocator, .energy),
        .snapshot_count = snapshotCount(config.steps, interval),
    };
}

pub fn runCavity(
    comptime dim: u8,
    allocator: std.mem.Allocator,
    config: CavityConfig(dim),
    writer: *std.Io.Writer,
) !RunResult(dim) {
    return switch (dim) {
        2 => try runCavity2D(allocator, config, writer),
        3 => try runCavity3D(allocator, config, writer),
        else => unreachable,
    };
}

fn runCavity2D(
    allocator: std.mem.Allocator,
    config: CavityConfig(2),
    writer: *std.Io.Writer,
) !RunResult(2) {
    const SystemType = System(2);
    const dt = config.dt();
    var mesh = try runtime.Mesh2D.cartesian(
        allocator,
        .{ config.grid, config.grid },
        .{ config.domain, config.domain },
    );
    defer mesh.deinit(allocator);

    try writer.writeAll("\n  ── TE₁₀ Cavity Resonance ──────────────────\n\n");
    try writer.print(
        "  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n\n",
        .{
            config.domain,
            config.grid,
            config.grid,
            2 * @as(u64, config.grid) * config.grid,
            config.spacing(),
            dt,
            config.courant,
        },
    );

    var system = try SystemType.cavity(allocator, &mesh, .{
        .domain_length = config.domain,
        .time_step = dt,
        .boundary = .pec,
    });
    defer system.deinit(allocator);

    const interval = @max(@as(u32, 1), common.framesToInterval(config.steps, config.frames));
    var evolution = try flux.evolution.Evolution(*SystemType, flux.integrators.Leapfrog).config()
        .dt(dt)
        .steps(config.steps)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(*SystemType)
                .field(.electric)
                .field(.magnetic)
                .measurement(.energy)
                .directory(config.output_dir)
                .baseName("cavity")
                .everySteps(interval),
        )
        .init(allocator, &system);
    defer evolution.deinit();

    const run_result = try evolution.run();
    return .{
        .elapsed_s = run_result.elapsed_s,
        .energy_final = try system.measurement(allocator, .energy),
        .snapshot_count = snapshotCount(config.steps, interval),
    };
}

fn runCavity3D(
    allocator: std.mem.Allocator,
    config: CavityConfig(3),
    writer: *std.Io.Writer,
) !RunResult(3) {
    const SystemType = System(3);
    var mesh = try runtime.Mesh3D.cartesian(
        allocator,
        .{ config.nx, config.ny, config.nz },
        .{ config.width, config.height, config.depth },
    );
    defer mesh.deinit(allocator);

    var system = try SystemType.cavity(allocator, &mesh, .{
        .width = config.width,
        .height = config.height,
        .time_step = config.dt,
        .boundary = .pec,
    });
    defer system.deinit(allocator);

    const omega = reference.tm110AngularFrequency(config.width, config.height);
    try writer.print(
        "\n  ── TM₁₁₀ Cavity Resonance (3D) ─────────────\n\n  domain    [0, {d:.2}] × [0, {d:.2}] × [0, {d:.2}]\n  grid      {d}×{d}×{d} ({d} tetrahedra)\n  spacing   h_min = {d:.6}\n  timestep  dt = {d:.6}\n  mode      TM₁₁₀  (ω = {d:.6})\n\n",
        .{
            config.width,
            config.height,
            config.depth,
            config.nx,
            config.ny,
            config.nz,
            mesh.num_tets(),
            config.gridSpacingMin(),
            config.dt,
            omega,
        },
    );

    const run_result = if (config.output_dir != null and config.output_interval > 0) blk: {
        var evolution = try flux.evolution.Evolution(*SystemType, flux.integrators.Leapfrog).config()
            .dt(config.dt)
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .listen(
                flux.listeners.Snapshots(*SystemType)
                    .field(.electric)
                    .field(.magnetic)
                    .measurement(.energy)
                    .directory(config.output_dir.?)
                    .baseName("maxwell_3d")
                    .everySteps(config.output_interval),
            )
            .init(allocator, &system);
        defer evolution.deinit();
        break :blk try evolution.run();
    } else blk: {
        var evolution = try flux.evolution.Evolution(*SystemType, flux.integrators.Leapfrog).config()
            .dt(config.dt)
            .steps(config.steps)
            .listen(flux.listeners.Progress(writer))
            .init(allocator, &system);
        defer evolution.deinit();
        break :blk try evolution.run();
    };

    _ = try reference.divergenceNorm3D(allocator, &system.state);
    return .{
        .elapsed_s = run_result.elapsed_s,
        .snapshot_count = if (config.output_dir != null and config.output_interval > 0)
            snapshotCount(config.steps, config.output_interval)
        else
            0,
    };
}

fn snapshotCount(steps: u32, interval: u32) u32 {
    std.debug.assert(interval > 0);
    const trailing: u32 = if (steps % interval == 0) 0 else 1;
    return 1 + (steps / interval) + trailing;
}

test {
    _ = @import("tests.zig");
}
