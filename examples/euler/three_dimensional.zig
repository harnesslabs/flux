const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const poisson = flux.operators.poisson;
const bridges = flux.operators.bridges;
const dec_context_mod = flux.operators.dec.context;
const feec_context_mod = flux.operators.feec.context;
const observers = flux.operators.observers;
const wedge_product = flux.operators.wedge_product;
const evolution_mod = flux.evolution;

pub const Mesh = flux.topology.Mesh(3, 3);
pub const Velocity = flux.forms.Cochain(Mesh, 1, flux.forms.Primal);
pub const Vorticity = flux.forms.Cochain(Mesh, 2, flux.forms.Primal);

pub const ConfigImpl = struct {
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
};

pub const RunResultImpl = struct {
    elapsed_s: f64,
    helicity_initial: f64,
    helicity_final: f64,
    snapshot_count: u32,
};

pub const SystemImpl = struct {
    mesh: *const Mesh,
    dec_operators: *dec_context_mod.OperatorContext(Mesh),
    feec_operators: *feec_context_mod.OperatorContext(Mesh),
    velocity: Velocity,
    vorticity: Vorticity,
    boundary_velocity: []f64,
    velocity_forcing: []f64,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh) !SystemImpl {
        var velocity = try Velocity.init(allocator, mesh);
        errdefer velocity.deinit(allocator);

        var vorticity = try Vorticity.init(allocator, mesh);
        errdefer vorticity.deinit(allocator);

        const dec_operators = try dec_context_mod.OperatorContext(Mesh).init(allocator, mesh);
        errdefer dec_operators.deinit();
        const feec_operators = try feec_context_mod.OperatorContext(Mesh).init(allocator, mesh);
        errdefer feec_operators.deinit();
        _ = try feec_operators.codifferential(2);
        _ = try feec_operators.codifferential(3);
        _ = try dec_operators.exteriorDerivative(flux.forms.Primal, 1);
        _ = try feec_operators.laplacian(1);

        const boundary_velocity = try allocator.alloc(f64, mesh.num_edges());
        errdefer allocator.free(boundary_velocity);
        @memset(boundary_velocity, 0.0);

        const velocity_forcing = try allocator.alloc(f64, mesh.num_edges());
        errdefer allocator.free(velocity_forcing);
        @memset(velocity_forcing, 0.0);

        return .{
            .mesh = mesh,
            .dec_operators = dec_operators,
            .feec_operators = feec_operators,
            .velocity = velocity,
            .vorticity = vorticity,
            .boundary_velocity = boundary_velocity,
            .velocity_forcing = velocity_forcing,
        };
    }

    pub fn deinit(self: *SystemImpl, allocator: std.mem.Allocator) void {
        allocator.free(self.velocity_forcing);
        allocator.free(self.boundary_velocity);
        self.vorticity.deinit(allocator);
        self.velocity.deinit(allocator);
        self.feec_operators.deinit();
        self.dec_operators.deinit();
    }
};

pub fn runImpl(allocator: std.mem.Allocator, config: ConfigImpl, writer: anytype) !RunResultImpl {
    var mesh = try Mesh.uniform_tetrahedral_grid(
        allocator,
        config.nx,
        config.ny,
        config.nz,
        config.width,
        config.height,
        config.depth,
    );
    defer mesh.deinit(allocator);

    var system = try SystemImpl.init(allocator, &mesh);
    defer system.deinit(allocator);
    try seedReferenceMode(allocator, &system);

    const helicity_initial = try conservedQuantity(allocator, &system);
    var evolution = try evolution_mod.Evolution(*SystemImpl, Euler3DMethod).config()
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
            .output_base_name = "euler_3d",
            .plan_override = if (config.output_dir != null and config.output_interval > 0)
                common.Plan.fromInterval(config.steps, config.output_interval, .{})
            else
                common.Plan.disabled(),
            .progress_writer = writer,
        },
        Renderer{ .state = &system },
    );

    const helicity_final = try conservedQuantity(allocator, &system);
    try writer.print(
        "euler_3d: grid={d}x{d}x{d} steps={d} dt={d:.6} helicity={d:.12} -> {d:.12}\n",
        .{ config.nx, config.ny, config.nz, config.steps, config.dt, helicity_initial, helicity_final },
    );

    return .{
        .elapsed_s = loop_result.elapsed_s,
        .helicity_initial = helicity_initial,
        .helicity_final = helicity_final,
        .snapshot_count = loop_result.snapshot_count,
    };
}

pub fn stepImpl(allocator: std.mem.Allocator, state: *SystemImpl, dt: f64) !void {
    _ = dt;
    try recoverVelocityFromVorticity(allocator, state);

    var vorticity = try (try state.dec_operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);

    const interpolate_velocity = bridges.WhitneyInterpolation(Mesh, 1).init(state.mesh);
    const interpolate_vorticity = bridges.WhitneyInterpolation(Mesh, 2).init(state.mesh);
    const advection_density = try wedge_product.wedge(
        allocator,
        interpolate_velocity.apply(&state.velocity),
        interpolate_vorticity.apply(&state.vorticity),
    );
    const project_density = bridges.DeRhamProjection(@TypeOf(advection_density).SpaceT).init(@TypeOf(advection_density).SpaceT.init(state.mesh));
    var projected_density = try project_density.apply(allocator, advection_density);
    defer projected_density.deinit(allocator);

    var transport = try (try state.feec_operators.codifferential(3)).apply(allocator, projected_density);
    defer transport.deinit(allocator);
}

pub fn seedReferenceMode(allocator: std.mem.Allocator, state: *SystemImpl) !void {
    const edge_vertices = state.mesh.simplices(1).items(.vertices);
    const coords = state.mesh.vertices.slice().items(.coords);
    for (state.velocity.values, edge_vertices) |*value, edge| {
        const p0 = coords[edge[0]];
        const p1 = coords[edge[1]];
        const midpoint = [3]f64{
            0.5 * (p0[0] + p1[0]),
            0.5 * (p0[1] + p1[1]),
            0.5 * (p0[2] + p1[2]),
        };
        const tangent = [3]f64{
            p1[0] - p0[0],
            p1[1] - p0[1],
            p1[2] - p0[2],
        };
        const envelope = midpoint[0] * (1.0 - midpoint[0]) *
            midpoint[1] * (1.0 - midpoint[1]) *
            midpoint[2] * (1.0 - midpoint[2]);
        const field = [3]f64{
            envelope * (@sin(std.math.pi * midpoint[2]) + 0.25 * @sin(2.0 * std.math.pi * midpoint[1])),
            envelope * (@sin(std.math.pi * midpoint[0]) + 0.20 * @sin(2.0 * std.math.pi * midpoint[2])),
            envelope * (@sin(std.math.pi * midpoint[1]) + 0.15 * @sin(2.0 * std.math.pi * midpoint[0])),
        };
        value.* = field[0] * tangent[0] + field[1] * tangent[1] + field[2] * tangent[2];
    }

    var vorticity = try (try state.dec_operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);
    @memcpy(state.boundary_velocity, state.velocity.values);

    var forcing = try (try state.feec_operators.laplacian(1)).apply(allocator, state.velocity);
    defer forcing.deinit(allocator);
    @memcpy(state.velocity_forcing, forcing.values);
}

pub fn conservedQuantity(allocator: std.mem.Allocator, state: *const SystemImpl) !f64 {
    const Helicity = observers.HelicityObserver(SystemImpl, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    return observer.evaluate(allocator, state, 0);
}

const Euler3DMethod = struct {
    pub fn advance(allocator: std.mem.Allocator, state: *SystemImpl, dt: f64) !void {
        try stepImpl(allocator, state, dt);
    }
};

fn selectVelocity(state: *const SystemImpl) *const Velocity {
    return &state.velocity;
}

fn recoverVelocityFromVorticity(allocator: std.mem.Allocator, state: *SystemImpl) !void {
    var solve = try poisson.solve_one_form_dirichlet(
        Mesh,
        allocator,
        state.feec_operators,
        state.velocity_forcing,
        state.boundary_velocity,
        .{
            .tolerance_relative = 1e-12,
            .iteration_limit = 4000,
        },
    );
    defer solve.deinit(allocator);

    @memcpy(state.velocity.values, solve.solution);
}

const Renderer = struct {
    state: *const SystemImpl,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try common.viz.writeProjectedTetFields(
            2,
            allocator,
            writer,
            self.state.mesh,
            .{
                .{ .name = "velocity_intensity", .kind = .edge_abs_mean, .values = self.state.velocity.values },
                .{ .name = "vorticity_flux", .kind = .face_abs_mean, .values = self.state.vorticity.values },
            },
        );
    }
};
