const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const common = @import("examples_common");
const core_mod = @import("maxwell_core");

const cochain = flux.forms;
const topology = flux.topology;
const operator_context_mod = flux.operators.context;

pub const Mesh3D = topology.Mesh(3, 3);

fn hooks(comptime MeshType: type) type {
    return struct {
        pub const apply_boundary_in_leapfrog = true;

        pub fn primeOperators(operators: *operator_context_mod.OperatorContext(MeshType)) !void {
            _ = try operators.exteriorDerivative(cochain.Primal, 1);
            _ = try operators.exteriorDerivative(cochain.Primal, 2);
            _ = try operators.exteriorDerivative(cochain.Dual, 1);
            _ = try operators.hodgeStar(2);
            _ = try operators.hodgeStarInverse(1);
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var star_b = try (try state.operators.hodgeStar(2)).apply(allocator, state.B);
            defer star_b.deinit(allocator);

            var derivative = try (try state.operators.exteriorDerivative(cochain.Dual, 1)).apply(allocator, star_b);
            defer derivative.deinit(allocator);

            var curl_b = try (try state.operators.hodgeStarInverse(1)).apply(allocator, derivative);
            defer curl_b.deinit(allocator);

            for (state.E.values, curl_b.values, state.J.values) |*electric, curl_value, current| {
                electric.* += dt * (curl_value - current);
            }
        }
    };
}

pub const Config = struct {
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

    fn gridSpacingMin(self: Config) f64 {
        const dx = self.width / @as(f64, @floatFromInt(self.nx));
        const dy = self.height / @as(f64, @floatFromInt(self.ny));
        const dz = self.depth / @as(f64, @floatFromInt(self.nz));
        return @min(dx, @min(dy, dz));
    }
};

pub fn State(comptime MeshType: type) type {
    return core_mod.Core(MeshType, hooks(MeshType)).State;
}

pub const MaxwellState3D = State(Mesh3D);

fn tm110WaveNumbers(width: f64, height: f64) struct { kx: f64, ky: f64 } {
    return .{
        .kx = std.math.pi / width,
        .ky = std.math.pi / height,
    };
}

fn tm110AngularFrequency(width: f64, height: f64) f64 {
    const wave = tm110WaveNumbers(width, height);
    return std.math.sqrt(wave.kx * wave.kx + wave.ky * wave.ky);
}

pub fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try core_mod.Core(@TypeOf(state.mesh.*), hooks(@TypeOf(state.mesh.*))).faradayStep(allocator, state, dt);
}

pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try core_mod.Core(@TypeOf(state.mesh.*), hooks(@TypeOf(state.mesh.*))).ampereStep(allocator, state, dt);
}

pub fn applyPecBoundary(_: std.mem.Allocator, state: anytype) !void {
    core_mod.Core(@TypeOf(state.mesh.*), hooks(@TypeOf(state.mesh.*))).applyPecBoundary(state);
}

pub fn leapfrogStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try core_mod.Core(@TypeOf(state.mesh.*), hooks(@TypeOf(state.mesh.*))).leapfrogStep(allocator, state, dt);
}

const SimResult = struct {
    elapsed_s: f64,
    snapshot_count: u32,
};

pub fn runSimulation(
    allocator: std.mem.Allocator,
    state: *MaxwellState3D,
    config: Config,
    writer: anytype,
) !SimResult {
    const plan: common.Plan = if (config.output_dir != null and config.output_interval > 0)
        common.Plan.fromInterval(config.steps, config.output_interval, .{})
    else
        common.Plan.disabled();

    var series = try common.Series.init(
        allocator,
        config.output_dir orelse "",
        "maxwell_3d",
        plan,
    );
    defer series.deinit();

    var progress = common.Progress(@TypeOf(writer)).init(writer, config.steps);

    var step_index: u32 = 0;
    while (step_index < config.steps) : (step_index += 1) {
        try leapfrogStep(allocator, state, config.dt);

        if (series.dueAt(step_index + 1, config.steps)) {
            const t = @as(f64, @floatFromInt(step_index + 1)) * config.dt;
            try series.capture(t, Maxwell3DRenderer{ .state = state });
        }

        progress.update(step_index + 1);
    }
    progress.finish();

    try series.finalize();

    return .{
        .elapsed_s = progress.elapsed(),
        .snapshot_count = series.count,
    };
}

const Maxwell3DRenderer = struct {
    state: *const MaxwellState3D,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try writeSnapshot(allocator, writer, self.state);
    }
};

/// Project the analytical TM₁₁₀ cavity electric field onto mesh edges.
///
/// TM₁₁₀ on a PEC box [0,a]×[0,b]×[0,c] with c = 1 has:
///   E = (0, 0, sin(πx/a) sin(πy/b) sin(ωt))
///   B = ( (k_y/ω) sin(k_x x) cos(k_y y) cos(ωt),
///        -(k_x/ω) cos(k_x x) sin(k_y y) cos(ωt),
///         0 )
/// with ω² = k_x² + k_y². This mode is uniform in z and still satisfies
/// PEC because E is normal to the z-walls and vanishes on x/y walls.
pub fn project_tm110_e(
    mesh: *const Mesh3D,
    values: []f64,
    t: f64,
    width: f64,
    height: f64,
) void {
    const wave = tm110WaveNumbers(width, height);
    const omega = tm110AngularFrequency(width, height);
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);

    for (values, edge_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const my = 0.5 * (p0[1] + p1[1]);
        const dz = p1[2] - p0[2];

        value.* = @sin(wave.kx * mx) * @sin(wave.ky * my) * @sin(omega * t) * dz;
    }
}

pub fn project_tm110_potential(
    mesh: *const Mesh3D,
    values: []f64,
    t: f64,
    width: f64,
    height: f64,
) void {
    const wave = tm110WaveNumbers(width, height);
    const omega = tm110AngularFrequency(width, height);
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);

    for (values, edge_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const my = 0.5 * (p0[1] + p1[1]);
        const dz = p1[2] - p0[2];

        value.* = (1.0 / omega) * @sin(wave.kx * mx) * @sin(wave.ky * my) * @cos(omega * t) * dz;
    }
}

/// Project the analytical TM₁₁₀ magnetic field onto mesh faces.
///
/// Each primal 2-form entry is a flux integral approximated by the face
/// centroid sample times the oriented face area vector.
pub fn project_tm110_b(
    mesh: *const Mesh3D,
    values: []f64,
    t: f64,
    width: f64,
    height: f64,
) void {
    const wave = tm110WaveNumbers(width, height);
    const omega = tm110AngularFrequency(width, height);
    const face_verts = mesh.simplices(2).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);

    for (values, face_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const p2 = coords[verts[2]];

        const centroid_x = (p0[0] + p1[0] + p2[0]) / 3.0;
        const centroid_y = (p0[1] + p1[1] + p2[1]) / 3.0;

        const bx = (wave.ky / omega) * @sin(wave.kx * centroid_x) * @cos(wave.ky * centroid_y) * @cos(omega * t);
        const by = -(wave.kx / omega) * @cos(wave.kx * centroid_x) * @sin(wave.ky * centroid_y) * @cos(omega * t);

        const edge_a = [3]f64{ p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        const edge_b = [3]f64{ p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        const area_vector = [3]f64{
            0.5 * (edge_a[1] * edge_b[2] - edge_a[2] * edge_b[1]),
            0.5 * (edge_a[2] * edge_b[0] - edge_a[0] * edge_b[2]),
            0.5 * (edge_a[0] * edge_b[1] - edge_a[1] * edge_b[0]),
        };

        value.* = bx * area_vector[0] + by * area_vector[1];
    }
}

pub fn seedTm110Mode(
    allocator: std.mem.Allocator,
    state: *MaxwellState3D,
    dt: f64,
    width: f64,
    height: f64,
) !void {
    @memset(state.E.values, 0.0);

    var potential = try MaxwellState3D.OneForm.init(allocator, state.mesh);
    defer potential.deinit(allocator);
    project_tm110_potential(state.mesh, potential.values, -dt / 2.0, width, height);

    var exact_flux = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, potential);
    defer exact_flux.deinit(allocator);
    @memcpy(state.B.values, exact_flux.values);
}

pub fn writeSnapshot(
    allocator: std.mem.Allocator,
    writer: anytype,
    state: *const MaxwellState3D,
) !void {
    try common.viz.writeProjectedTetFields(
        2,
        allocator,
        writer,
        state.mesh,
        .{
            .{ .name = "E_intensity", .kind = .edge_abs_mean, .values = state.E.values },
            .{ .name = "B_flux", .kind = .face_abs_mean, .values = state.B.values },
        },
    );
}

pub fn runDriver(allocator: std.mem.Allocator, config: Config) !void {
    var mesh = try makeCavityMesh(allocator, config);
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedTm110Mode(allocator, &state, config.dt, config.width, config.height);
    const omega = tm110AngularFrequency(config.width, config.height);

    var stderr_buffer: [1024]u8 = undefined;
    var stderr_writer = std.fs.File.stderr().writer(&stderr_buffer);
    const writer = &stderr_writer.interface;

    try writer.writeAll("\n  ── TM₁₁₀ Cavity Resonance (3D) ─────────────\n\n");
    try writer.print("  domain    [0, {d:.2}] × [0, {d:.2}] × [0, {d:.2}]\n", .{
        config.width, config.height, config.depth,
    });
    try writer.print("  grid      {d}×{d}×{d} ({d} tetrahedra)\n", .{
        config.nx, config.ny, config.nz, mesh.num_tets(),
    });
    try writer.print("  spacing   h_min = {d:.6}\n", .{config.gridSpacingMin()});
    try writer.print("  timestep  dt = {d:.6}\n", .{config.dt});
    try writer.print("  mode      TM₁₁₀  (ω = {d:.6})\n", .{omega});
    try writer.print("  mesh      {d} vertices  {d} edges  {d} faces  {d} tets\n\n", .{
        mesh.num_vertices(), mesh.num_edges(), mesh.num_faces(), mesh.num_tets(),
    });

    const result = try runSimulation(allocator, &state, config, writer);

    const divergence = try divergenceNorm(allocator, &state);
    const steps_per_sec = @as(f64, @floatFromInt(config.steps)) / result.elapsed_s;
    var duration_buf: [16]u8 = undefined;

    try writer.writeAll("\n  ── Results ─────────────────────────────────\n\n");
    try writer.print("  steps    {d}\n", .{state.timestep});
    try writer.print("  omega    {d:.6}\n", .{omega});
    try writer.print("  ||dB||₂  {e}\n", .{divergence});
    try writer.print("  elapsed  {s} ({d:.0} steps/s)\n", .{
        common.formatDuration(&duration_buf, result.elapsed_s),
        steps_per_sec,
    });
    if (result.snapshot_count > 0 and config.output_dir != null) {
        try writer.print("  output   {d} frames → {s}/\n", .{
            result.snapshot_count,
            config.output_dir.?,
        });
        try writer.print("\n  ▸ uv run tools/visualize.py {s} --field B_flux --output {s}/animation.png\n\n", .{
            config.output_dir.?,
            config.output_dir.?,
        });
    } else {
        try writer.writeAll("\n");
    }
}

pub fn makeCavityMesh(allocator: std.mem.Allocator, config: Config) !Mesh3D {
    return Mesh3D.uniform_tetrahedral_grid(
        allocator,
        config.nx,
        config.ny,
        config.nz,
        config.width,
        config.height,
        config.depth,
    );
}

fn seedClosedMagneticField(allocator: std.mem.Allocator, state: *MaxwellState3D) !void {
    var potential = try MaxwellState3D.OneForm.init(allocator, state.mesh);
    defer potential.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0x92_3D_B_00);
    for (potential.values) |*value| {
        value.* = rng.random().float(f64) * 0.2 - 0.1;
    }

    var exact_flux = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, potential);
    defer exact_flux.deinit(allocator);

    @memcpy(state.B.values, exact_flux.values);
}

pub fn divergenceNorm(allocator: std.mem.Allocator, state: *const MaxwellState3D) !f64 {
    var divergence = try (try state.operators.exteriorDerivative(cochain.Primal, 2)).apply(allocator, state.B);
    defer divergence.deinit(allocator);
    return std.math.sqrt(divergence.norm_squared());
}

fn compute_tm110_eigenvalue(
    allocator: std.mem.Allocator,
    nx: u32,
    ny: u32,
    nz: u32,
    width: f64,
    height: f64,
    depth: f64,
) !f64 {
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, nx, ny, nz, width, height, depth);
    defer mesh.deinit(allocator);

    const operator_context = try operator_context_mod.OperatorContext(Mesh3D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.exteriorDerivative(cochain.Primal, 1);
    _ = try operator_context.hodgeStar(1);
    _ = try operator_context.hodgeStar(2);

    var E = try MaxwellState3D.OneForm.init(allocator, &mesh);
    defer E.deinit(allocator);
    const omega = tm110AngularFrequency(width, height);
    project_tm110_e(&mesh, E.values, std.math.pi / (2.0 * omega), width, height);

    var dE = try (try operator_context.exteriorDerivative(cochain.Primal, 1)).apply(allocator, E);
    defer dE.deinit(allocator);

    var star_dE = try (try operator_context.hodgeStar(2)).apply(allocator, dE);
    defer star_dE.deinit(allocator);
    var numerator: f64 = 0.0;
    for (dE.values, star_dE.values) |lhs, rhs| {
        numerator += lhs * rhs;
    }

    var star_E = try (try operator_context.hodgeStar(1)).apply(allocator, E);
    defer star_E.deinit(allocator);
    var denominator: f64 = 0.0;
    for (E.values, star_E.values) |lhs, rhs| {
        denominator += lhs * rhs;
    }

    return numerator / denominator;
}

test "MaxwellState3D initializes field spaces for tetrahedral meshes" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 1, .ny = 1, .nz = 1 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(@as(u64, 0), state.timestep);
    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.E.values.len)));
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(state.B.values.len)));
    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.J.values.len)));

    comptime {
        try testing.expectEqual(1, MaxwellState3D.OneForm.degree);
        try testing.expectEqual(2, MaxwellState3D.TwoForm.degree);
    }
}

test "3D Maxwell preserves d₂B = 0 over 1000 source-free cavity steps" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 2, .ny = 2, .nz = 2, .dt = 0.0025 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedTm110Mode(allocator, &state, 0.0025, 1.0, 1.0);

    for (0..1000) |_| {
        try leapfrogStep(allocator, &state, 0.0025);

        const norm = try divergenceNorm(allocator, &state);
        try testing.expectApproxEqAbs(@as(f64, 0.0), norm, 1e-12);
    }
}

test "3D Maxwell PEC zeroes boundary edge circulation after a step" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 2, .ny = 1, .nz = 1, .dt = 0.0025 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    for (state.E.values, 0..) |*value, idx| {
        value.* = @as(f64, @floatFromInt(idx + 1));
    }

    try applyPecBoundary(allocator, &state);

    for (mesh.boundary_edges) |edge_idx| {
        try testing.expectEqual(@as(f64, 0.0), state.E.values[edge_idx]);
    }
}

test "3D Maxwell snapshot exports projected E and B data on tetrahedral cells" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 1, .ny = 1, .nz = 1 });
    defer mesh.deinit(allocator);

    var state = try MaxwellState3D.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedTm110Mode(allocator, &state, 0.0025, 1.0, 1.0);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try writeSnapshot(allocator, output.writer(allocator), &state);

    const xml = output.items;
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"E_intensity\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"B_flux\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "NumberOfCells=\"6\"") != null);
}

test "project_tm110_e vanishes on non-z edges and x/y boundary z-edges" {
    const allocator = testing.allocator;
    var mesh = try makeCavityMesh(allocator, .{ .nx = 2, .ny = 2, .nz = 2 });
    defer mesh.deinit(allocator);

    const values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(values);
    const omega = tm110AngularFrequency(1.0, 1.0);
    project_tm110_e(&mesh, values, std.math.pi / (2.0 * omega), 1.0, 1.0);

    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_verts) |value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const dz = p1[2] - p0[2];
        const mx = 0.5 * (p0[0] + p1[0]);
        const my = 0.5 * (p0[1] + p1[1]);
        if (@abs(dz) < 1e-15) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-15);
        }
        if (@abs(mx) < 1e-12 or @abs(mx - 1.0) < 1e-12 or @abs(my) < 1e-12 or @abs(my - 1.0) < 1e-12) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-14);
        }
    }
}

test "TM₁₁₀ cavity: eigenvalue error decreases under 3D refinement" {
    const allocator = testing.allocator;
    const analytical = 2.0 * std.math.pi * std.math.pi;

    const coarse = try compute_tm110_eigenvalue(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    const fine = try compute_tm110_eigenvalue(allocator, 4, 4, 4, 1.0, 1.0, 1.0);

    const error_coarse = @abs(coarse - analytical) / analytical;
    const error_fine = @abs(fine - analytical) / analytical;

    std.debug.print("\nTM110 ω²: 2^3={d:.6}, 4^3={d:.6}, analytical={d:.6}\n", .{
        coarse, fine, analytical,
    });
    std.debug.print("TM110 error: 2^3={d:.6}, 4^3={d:.6}, ratio={d:.2}\n", .{
        error_coarse, error_fine, error_coarse / error_fine,
    });

    try testing.expect(error_fine < error_coarse);
    try testing.expect(error_coarse / error_fine >= 1.5);
}
