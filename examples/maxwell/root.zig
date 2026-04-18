const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const runtime = @import("runtime.zig");
const reference = @import("reference.zig");

const Demo = enum { dipole, cavity };

const BoundaryCondition = enum {
    pec,
};

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

fn DipoleOptions(comptime dim: u8) type {
    return struct {
        center: [dim]f64,
        frequency: f64,
        amplitude: f64,
        boundary: BoundaryCondition = .pec,
    };
}

fn CavityOptions2D() type {
    return struct {
        domain_length: f64,
        time_step: f64,
        boundary: BoundaryCondition = .pec,
    };
}

fn CavityOptions3D() type {
    return struct {
        width: f64,
        height: f64,
        time_step: f64,
        boundary: BoundaryCondition = .pec,
    };
}

fn Maxwell2D(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const Field = enum { electric, magnetic, current };
        pub const Measurement = enum { energy };

        state: runtime.StateForMesh2D(MeshType),
        source: ?runtime.PointDipole(MeshType) = null,
        boundary: BoundaryCondition = .pec,

        pub const Leapfrog = struct {
            pub fn first(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                const time = @as(f64, @floatFromInt(system.state.timestep)) * dt;
                if (system.source) |source| {
                    source.apply(&system.state.J, time);
                } else {
                    @memset(system.state.J.values, 0.0);
                }
                try runtime.faraday_step(allocator, &system.state, dt);
            }

            pub fn second(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try runtime.ampere_step(allocator, &system.state, dt);
                system.state.timestep += 1;
            }

            pub fn applyBoundary(system: *Self) void {
                switch (system.boundary) {
                    .pec => runtime.apply_pec_boundary(&system.state),
                }
            }
        };

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = try runtime.StateForMesh2D(MeshType).init(allocator, mesh),
            };
        }

        pub fn dipole(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: DipoleOptions(2),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            system.source = runtime.PointDipole(MeshType).init(mesh, options.frequency, options.amplitude, options.center);
            system.boundary = options.boundary;
            return system;
        }

        pub fn cavity(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: CavityOptions2D(),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            reference.project_te10_b(mesh, system.state.B.values, -options.time_step / 2.0, options.domain_length);
            system.boundary = options.boundary;
            system.source = null;
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return switch (which) {
                .energy => try runtime.electromagnetic_energy(allocator, &self.state),
            };
        }

        pub fn writeFields(self: *const Self, allocator: std.mem.Allocator, writer: anytype, fields: []const Field) !void {
            var cell_data: [3]flux.io.DataArraySlice = undefined;
            var projected: [3]?[]f64 = .{ null, null, null };
            defer for (projected) |maybe_values| {
                if (maybe_values) |values| allocator.free(values);
            };

            var count: usize = 0;
            for (fields) |field| {
                switch (field) {
                    .electric => {
                        projected[count] = try flux.io.project_edges_to_faces(allocator, self.state.mesh.*, self.state.E.values);
                        cell_data[count] = .{ .name = "electric", .values = projected[count].? };
                    },
                    .magnetic => {
                        cell_data[count] = .{ .name = "magnetic", .values = self.state.B.values };
                    },
                    .current => {
                        projected[count] = try flux.io.project_edges_to_faces(allocator, self.state.mesh.*, self.state.J.values);
                        cell_data[count] = .{ .name = "current", .values = projected[count].? };
                    },
                }
                count += 1;
            }

            try flux.io.write(writer, self.state.mesh.*, &.{}, cell_data[0..count]);
        }
    };
}

fn Maxwell3D(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const Field = enum { electric, magnetic, current };
        pub const Measurement = enum { energy };

        state: runtime.StateForMesh3D(MeshType),
        boundary: BoundaryCondition = .pec,

        pub const Leapfrog = struct {
            pub fn first(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try runtime.faraday_step(allocator, &system.state, dt);
            }

            pub fn second(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try runtime.ampere_step_3d(allocator, &system.state, dt);
                system.state.timestep += 1;
            }

            pub fn applyBoundary(system: *Self) void {
                switch (system.boundary) {
                    .pec => runtime.apply_pec_boundary(&system.state),
                }
            }
        };

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = try runtime.StateForMesh3D(MeshType).init(allocator, mesh),
            };
        }

        pub fn cavity(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: CavityOptions3D(),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            try reference.seedTm110Mode(allocator, &system.state, options.time_step, options.width, options.height);
            system.boundary = options.boundary;
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return switch (which) {
                .energy => try runtime.electromagnetic_energy(allocator, &self.state),
            };
        }

        pub fn writeFields(self: *const Self, allocator: std.mem.Allocator, writer: anytype, fields: []const Field) !void {
            var projected_fields: [3]common.viz.TetProjectionField = undefined;
            var count: usize = 0;
            for (fields) |field| {
                projected_fields[count] = switch (field) {
                    .electric => .{ .name = "electric", .kind = .edge_abs_mean, .values = self.state.E.values },
                    .magnetic => .{ .name = "magnetic", .kind = .face_abs_mean, .values = self.state.B.values },
                    .current => .{ .name = "current", .kind = .edge_abs_mean, .values = self.state.J.values },
                };
                count += 1;
            }

            switch (count) {
                0 => try flux.io.write(writer, self.state.mesh.*, &.{}, &.{}),
                1 => try common.viz.writeProjectedTetFields(1, allocator, writer, self.state.mesh, .{projected_fields[0]}),
                2 => try common.viz.writeProjectedTetFields(2, allocator, writer, self.state.mesh, .{
                    projected_fields[0],
                    projected_fields[1],
                }),
                3 => try common.viz.writeProjectedTetFields(3, allocator, writer, self.state.mesh, .{
                    projected_fields[0],
                    projected_fields[1],
                    projected_fields[2],
                }),
                else => unreachable,
            }
        }
    };
}

pub fn Maxwell(comptime dim: u8, comptime MeshType: type) type {
    return switch (dim) {
        2 => Maxwell2D(MeshType),
        3 => Maxwell3D(MeshType),
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

fn run2D(allocator: std.mem.Allocator, config: Config2D, writer: *std.Io.Writer) !RunResult2D {
    const SystemType = Maxwell(2, runtime.Mesh2D);
    const dt = config.dt();
    var mesh = try runtime.Mesh2D.plane(allocator, config.grid, config.grid, config.domain, config.domain);
    defer mesh.deinit(allocator);

    const base_name: []const u8 = switch (config.demo) {
        .dipole => "dipole",
        .cavity => "cavity",
    };

    var system = switch (config.demo) {
        .dipole => blk: {
            const frequency = config.sourceFrequency();
            try writer.writeAll("\n  ── Dipole Simulation ───────────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n  source    f = {d:.4} Hz, A = {d:.2}\n\n", .{
                config.domain,
                config.grid,
                config.grid,
                2 * @as(u64, config.grid) * config.grid,
                config.spacing(),
                dt,
                config.courant,
                frequency,
                config.amplitude,
            });
            break :blk try SystemType.dipole(allocator, &mesh, .{
                .center = .{ config.domain / 2.0, config.domain / 2.0 },
                .frequency = frequency,
                .amplitude = config.amplitude,
            });
        },
        .cavity => blk: {
            try writer.writeAll("\n  ── TE₁₀ Cavity Resonance ──────────────────\n\n");
            try writer.print("  domain    [0, {d:.2}]²\n  grid      {d}×{d} ({d} triangles)\n  spacing   h = {d:.6}\n  timestep  dt = {d:.6} (Courant {d:.2})\n\n", .{
                config.domain,
                config.grid,
                config.grid,
                2 * @as(u64, config.grid) * config.grid,
                config.spacing(),
                dt,
                config.courant,
            });
            break :blk try SystemType.cavity(allocator, &mesh, .{
                .domain_length = config.domain,
                .time_step = dt,
            });
        },
    };
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
                .baseName(base_name)
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

fn snapshotCount(steps: u32, interval: u32) u32 {
    std.debug.assert(interval > 0);
    const trailing: u32 = if (steps % interval == 0) 0 else 1;
    return 1 + (steps / interval) + trailing;
}

fn run3D(allocator: std.mem.Allocator, config: Config3D, writer: *std.Io.Writer) !RunResult3D {
    const SystemType = Maxwell(3, runtime.Mesh3D);
    var mesh = try makeCavityMesh(allocator, config);
    defer mesh.deinit(allocator);

    var system = try SystemType.cavity(allocator, &mesh, .{
        .width = config.width,
        .height = config.height,
        .time_step = config.dt,
    });
    defer system.deinit(allocator);

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
    return Maxwell(dim, Mesh(dim));
}

pub fn run(comptime dim: u8, allocator: std.mem.Allocator, config: Config(dim), writer: *std.Io.Writer) !RunResult(dim) {
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
    const Method = flux.integrators.Leapfrog(*System(dim));
    try Method.advance(allocator, system, dt);
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
