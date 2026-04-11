const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const runtime = @import("runtime.zig");
const reference = @import("reference.zig");

const flux_io = flux.io;

const Demo = enum { dipole, cavity };

const Config2D = struct {
    demo: Demo = .dipole,
    steps: u32 = 1000,
    grid: u32 = 32,
    domain: f64 = 1.0,
    courant: f64 = 0.1,
    frequency: f64 = 0.0,
    amplitude: f64 = 1.0,
    output_dir: []const u8 = "output",
    frames: u32 = 100,

    pub fn spacing(self: Config2D) f64 {
        return self.domain / @as(f64, @floatFromInt(self.grid));
    }

    pub fn dt(self: Config2D) f64 {
        return self.courant * self.spacing();
    }

    pub fn sourceFrequency(self: Config2D) f64 {
        if (self.frequency != 0.0) return self.frequency;
        return 1.0 / (2.0 * self.domain);
    }
};

const RunResult2D = struct {
    elapsed_s: f64,
    energy_final: f64,
    snapshot_count: u32,
};

const Maxwell2DRenderer = struct {
    state: *const runtime.MaxwellState2D,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try flux_io.write_fields(allocator, writer, self.state.mesh.*, self.state.E.values, self.state.B.values);
    }
};

fn simulate2D(allocator: std.mem.Allocator, state: *runtime.MaxwellState2D, source: ?runtime.PointDipole(runtime.Mesh2D), config: Config2D, base_name: []const u8, writer: anytype) !RunResult2D {
    const dt = config.dt();
    var series = try common.Series.init(allocator, config.output_dir, base_name, common.Plan.fromFrames(config.steps, config.frames, .{}));
    defer series.deinit();

    var progress = common.Progress.init(writer, config.steps);
    for (0..config.steps) |step_idx| {
        const t = @as(f64, @floatFromInt(step_idx)) * dt;
        if (source) |dipole| dipole.apply(&state.J, t);
        try runtime.leapfrog_step(allocator, state, dt);
        runtime.apply_pec_boundary(state);
        if (series.dueAt(@intCast(step_idx + 1), config.steps)) {
            try series.capture(t + dt, Maxwell2DRenderer{ .state = state });
        }
        progress.update(@intCast(step_idx + 1));
    }
    progress.finish();
    try series.finalize();
    return .{
        .elapsed_s = progress.elapsed(),
        .energy_final = try runtime.electromagnetic_energy(allocator, state),
        .snapshot_count = series.count,
    };
}

fn run2D(allocator: std.mem.Allocator, config: Config2D, writer: anytype) !RunResult2D {
    const h = config.spacing();
    const dt = config.dt();
    var mesh = try runtime.Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);
    var state = try runtime.MaxwellState2D.init(allocator, &mesh);
    defer state.deinit(allocator);

    switch (config.demo) {
        .dipole => {
            const frequency = config.sourceFrequency();
            try writer.writeAll("\n  ── Dipole Simulation ───────────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n  source    f = {d:.4} Hz, A = {d:.2}\n\n", .{
                config.domain, config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid, h, dt, config.courant, frequency, config.amplitude,
            });
            const center = [2]f64{ config.domain / 2.0, config.domain / 2.0 };
            const dipole = runtime.PointDipole(runtime.Mesh2D).init(&mesh, frequency, config.amplitude, center);
            return simulate2D(allocator, &state, dipole, config, "dipole", writer);
        },
        .cavity => {
            try writer.writeAll("\n  ── TE₁₀ Cavity Resonance ──────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n\n", .{
                config.domain, config.grid, config.grid, 2 * @as(u64, config.grid) * config.grid, h, dt, config.courant,
            });
            reference.project_te10_b(&mesh, state.B.values, -dt / 2.0, config.domain);
            return simulate2D(allocator, &state, null, config, "cavity", writer);
        },
    }
}

const Config3D = struct {
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

    fn gridSpacingMin(self: Config3D) f64 {
        return @min(self.width / @as(f64, @floatFromInt(self.nx)), @min(self.height / @as(f64, @floatFromInt(self.ny)), self.depth / @as(f64, @floatFromInt(self.nz))));
    }
};

const RunResult3D = struct {
    elapsed_s: f64,
    snapshot_count: u32,
};

const Maxwell3DRenderer = struct {
    state: *const runtime.MaxwellState3D,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try reference.writeSnapshot(allocator, writer, self.state);
    }
};

fn makeCavityMesh(allocator: std.mem.Allocator, config: Config3D) !runtime.Mesh3D {
    return runtime.Mesh3D.uniform_tetrahedral_grid(allocator, config.nx, config.ny, config.nz, config.width, config.height, config.depth);
}

fn run3D(allocator: std.mem.Allocator, config: Config3D, writer: anytype) !RunResult3D {
    var mesh = try makeCavityMesh(allocator, config);
    defer mesh.deinit(allocator);
    var state = try runtime.MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);
    try reference.seedTm110Mode(allocator, &state, config.dt, config.width, config.height);

    const omega = reference.tm110AngularFrequency(config.width, config.height);
    try writer.print("\n  ── TM₁₁₀ Cavity Resonance (3D) ─────────────\n\n  domain    [0, {d:.2}] × [0, {d:.2}] × [0, {d:.2}]\n  grid      {d}×{d}×{d} ({d} tetrahedra)\n  spacing   h_min = {d:.6}\n  timestep  dt = {d:.6}\n  mode      TM₁₁₀  (ω = {d:.6})\n\n", .{
        config.width, config.height, config.depth, config.nx, config.ny, config.nz, mesh.num_tets(), config.gridSpacingMin(), config.dt, omega,
    });

    const plan: common.Plan = if (config.output_dir != null and config.output_interval > 0) common.Plan.fromInterval(config.steps, config.output_interval, .{}) else common.Plan.disabled();
    var series = try common.Series.init(allocator, config.output_dir orelse "", "maxwell_3d", plan);
    defer series.deinit();
    var progress = common.Progress.init(writer, config.steps);
    for (0..config.steps) |step_idx| {
        try runtime.leapfrog_step_3d(allocator, &state, config.dt);
        if (series.dueAt(@intCast(step_idx + 1), config.steps)) {
            const t = @as(f64, @floatFromInt(step_idx + 1)) * config.dt;
            try series.capture(t, Maxwell3DRenderer{ .state = &state });
        }
        progress.update(@intCast(step_idx + 1));
    }
    progress.finish();
    try series.finalize();
    _ = try reference.divergenceNorm3D(allocator, &state);
    return .{ .elapsed_s = progress.elapsed(), .snapshot_count = series.count };
}

pub fn Mesh(comptime dim: u8) type {
    return switch (dim) {
        2 => runtime.Mesh2D,
        3 => runtime.Mesh3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn Config(comptime dim: u8) type {
    return switch (dim) {
        2 => Config2D,
        3 => Config3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn RunResult(comptime dim: u8) type {
    return switch (dim) {
        2 => RunResult2D,
        3 => RunResult3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn State(comptime dim: u8) type {
    return switch (dim) {
        2 => runtime.MaxwellState2D,
        3 => runtime.MaxwellState3D,
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

pub fn run(comptime dim: u8, allocator: std.mem.Allocator, config: Config(dim), writer: anytype) !RunResult(dim) {
    return switch (dim) {
        2 => try run2D(allocator, config, writer),
        3 => try run3D(allocator, config, writer),
        else => unreachable,
    };
}

pub fn makeMesh(comptime dim: u8, allocator: std.mem.Allocator, config: Config(dim)) !Mesh(dim) {
    return switch (dim) {
        2 => try runtime.Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain),
        3 => try makeCavityMesh(allocator, config),
        else => unreachable,
    };
}

pub fn step(comptime dim: u8, allocator: std.mem.Allocator, state: *State(dim), dt: f64) !void {
    switch (dim) {
        2 => {
            try runtime.leapfrog_step(allocator, state, dt);
            runtime.apply_pec_boundary(state);
        },
        3 => try runtime.leapfrog_step_3d(allocator, state, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(comptime dim: u8, allocator: std.mem.Allocator, state: *State(dim), dt: f64, width: f64, height: f64) !void {
    switch (dim) {
        2 => reference.project_te10_b(state.mesh, state.B.values, -dt / 2.0, width),
        3 => try reference.seedTm110Mode(allocator, state, dt, width, height),
        else => unreachable,
    }
}

pub fn divergenceNorm(allocator: std.mem.Allocator, state: *const State(3)) !f64 {
    return reference.divergenceNorm3D(allocator, state);
}

test {
    _ = @import("tests.zig");
}
