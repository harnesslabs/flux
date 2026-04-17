const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const runtime = @import("runtime.zig");
const reference = @import("reference.zig");

const flux_io = flux.io;
const evolution_mod = flux.evolution;

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

pub const System2D = struct {
    state: runtime.MaxwellState2D,
    source: ?runtime.PointDipole(runtime.Mesh2D) = null,

    pub fn init(allocator: std.mem.Allocator, mesh: *const runtime.Mesh2D) !System2D {
        return .{
            .state = try runtime.MaxwellState2D.init(allocator, mesh),
        };
    }

    pub fn deinit(self: *System2D, allocator: std.mem.Allocator) void {
        self.state.deinit(allocator);
    }
};

const Maxwell2DRenderer = struct {
    system: *const System2D,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try flux_io.write_fields(
            allocator,
            writer,
            self.system.state.mesh.*,
            self.system.state.E.values,
            self.system.state.B.values,
        );
    }
};

const DrivenLeapfrog2DMethod = struct {
    pub fn advance(
        allocator: std.mem.Allocator,
        system: *System2D,
        dt: f64,
    ) !void {
        const time = @as(f64, @floatFromInt(system.state.timestep)) * dt;
        if (system.source) |dipole| dipole.apply(&system.state.J, time);
        try runtime.leapfrog_step(allocator, &system.state, dt);
        runtime.apply_pec_boundary(&system.state);
    }
};

fn simulate2D(
    allocator: std.mem.Allocator,
    system: *System2D,
    config: Config2D,
    base_name: []const u8,
    writer: anytype,
) !RunResult2D {
    const dt = config.dt();
    var evolution = try evolution_mod.Evolution(*System2D, DrivenLeapfrog2DMethod).config()
        .dt(dt)
        .init(allocator, system);
    defer evolution.deinit();

    const loop_result = try common.runEvolutionLoop(
        allocator,
        &evolution,
        .{
            .steps = config.steps,
            .final_time = dt * @as(f64, @floatFromInt(config.steps)),
            .frames = config.frames,
            .output_dir = config.output_dir,
            .output_base_name = base_name,
            .progress_writer = writer,
        },
        Maxwell2DRenderer{ .system = system },
    );

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .energy_final = try runtime.electromagnetic_energy(allocator, &system.state),
        .snapshot_count = loop_result.snapshot_count,
    };
}

fn run2D(allocator: std.mem.Allocator, config: Config2D, writer: anytype) !RunResult2D {
    const h = config.spacing();
    const dt = config.dt();
    var mesh = try runtime.Mesh2D.plane(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    var system = try System2D.init(allocator, &mesh);
    defer system.deinit(allocator);

    switch (config.demo) {
        .dipole => {
            const frequency = config.sourceFrequency();
            try writer.writeAll("\n  ── Dipole Simulation ───────────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n  source    f = {d:.4} Hz, A = {d:.2}\n\n", .{
                config.domain,
                config.grid,
                config.grid,
                2 * @as(u64, config.grid) * config.grid,
                h,
                dt,
                config.courant,
                frequency,
                config.amplitude,
            });
            const center = [2]f64{ config.domain / 2.0, config.domain / 2.0 };
            system.source = runtime.PointDipole(runtime.Mesh2D).init(&mesh, frequency, config.amplitude, center);
            return simulate2D(allocator, &system, config, "dipole", writer);
        },
        .cavity => {
            try writer.writeAll("\n  ── TE₁₀ Cavity Resonance ──────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n\n", .{
                config.domain,
                config.grid,
                config.grid,
                2 * @as(u64, config.grid) * config.grid,
                h,
                dt,
                config.courant,
            });
            reference.project_te10_b(&mesh, system.state.B.values, -dt / 2.0, config.domain);
            system.source = null;
            return simulate2D(allocator, &system, config, "cavity", writer);
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
        return @min(
            self.width / @as(f64, @floatFromInt(self.nx)),
            @min(
                self.height / @as(f64, @floatFromInt(self.ny)),
                self.depth / @as(f64, @floatFromInt(self.nz)),
            ),
        );
    }
};

const RunResult3D = struct {
    elapsed_s: f64,
    snapshot_count: u32,
};

pub const System3D = struct {
    state: runtime.MaxwellState3D,

    pub fn init(allocator: std.mem.Allocator, mesh: *const runtime.Mesh3D) !System3D {
        return .{
            .state = try runtime.MaxwellState3D.init(allocator, mesh),
        };
    }

    pub fn deinit(self: *System3D, allocator: std.mem.Allocator) void {
        self.state.deinit(allocator);
    }
};

const Maxwell3DRenderer = struct {
    system: *const System3D,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try reference.writeSnapshot(allocator, writer, &self.system.state);
    }
};

const Leapfrog3DMethod = struct {
    pub fn advance(
        allocator: std.mem.Allocator,
        system: *System3D,
        dt: f64,
    ) !void {
        try runtime.leapfrog_step_3d(allocator, &system.state, dt);
    }
};

fn makeCavityMesh(allocator: std.mem.Allocator, config: Config3D) !runtime.Mesh3D {
    return runtime.Mesh3D.uniform_tetrahedral_grid(
        allocator,
        config.nx,
        config.ny,
        config.nz,
        config.width,
        config.height,
        config.depth,
    );
}

fn run3D(allocator: std.mem.Allocator, config: Config3D, writer: anytype) !RunResult3D {
    var mesh = try makeCavityMesh(allocator, config);
    defer mesh.deinit(allocator);

    var system = try System3D.init(allocator, &mesh);
    defer system.deinit(allocator);
    try reference.seedTm110Mode(allocator, &system.state, config.dt, config.width, config.height);

    const omega = reference.tm110AngularFrequency(config.width, config.height);
    try writer.print("\n  ── TM₁₁₀ Cavity Resonance (3D) ─────────────\n\n  domain    [0, {d:.2}] × [0, {d:.2}] × [0, {d:.2}]\n  grid      {d}×{d}×{d} ({d} tetrahedra)\n  spacing   h_min = {d:.6}\n  timestep  dt = {d:.6}\n  mode      TM₁₁₀  (ω = {d:.6})\n\n", .{
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
    });

    var evolution = try evolution_mod.Evolution(*System3D, Leapfrog3DMethod).config()
        .dt(config.dt)
        .init(allocator, &system);
    defer evolution.deinit();

    const loop_result = try common.runEvolutionLoop(
        allocator,
        &evolution,
        .{
            .steps = config.steps,
            .final_time = config.dt * @as(f64, @floatFromInt(config.steps)),
            .frames = 0,
            .output_dir = config.output_dir orelse "",
            .output_base_name = "maxwell_3d",
            .plan_override = if (config.output_dir != null and config.output_interval > 0)
                common.Plan.fromInterval(config.steps, config.output_interval, .{})
            else
                common.Plan.disabled(),
            .progress_writer = writer,
        },
        Maxwell3DRenderer{ .system = &system },
    );

    _ = try reference.divergenceNorm3D(allocator, &system.state);
    return .{ .elapsed_s = loop_result.elapsed_s, .snapshot_count = loop_result.snapshot_count };
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

pub fn System(comptime dim: u8) type {
    return switch (dim) {
        2 => System2D,
        3 => System3D,
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
        2 => try runtime.Mesh2D.plane(allocator, config.grid, config.grid, config.domain, config.domain),
        3 => try makeCavityMesh(allocator, config),
        else => unreachable,
    };
}

pub fn step(comptime dim: u8, allocator: std.mem.Allocator, system: *System(dim), dt: f64) !void {
    switch (dim) {
        2 => try DrivenLeapfrog2DMethod.advance(allocator, system, dt),
        3 => try Leapfrog3DMethod.advance(allocator, system, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(
    comptime dim: u8,
    allocator: std.mem.Allocator,
    system: *System(dim),
    dt: f64,
    width: f64,
    height: f64,
) !void {
    switch (dim) {
        2 => reference.project_te10_b(system.state.mesh, system.state.B.values, -dt / 2.0, width),
        3 => try reference.seedTm110Mode(allocator, &system.state, dt, width, height),
        else => unreachable,
    }
}

pub fn divergenceNorm(allocator: std.mem.Allocator, system: *const System(3)) !f64 {
    return reference.divergenceNorm3D(allocator, &system.state);
}

test {
    _ = @import("tests.zig");
}
