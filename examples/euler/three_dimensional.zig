const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const poisson = flux.operators.poisson;
const operator_context_mod = flux.operators.context;
const observers = flux.operators.observers;
const wedge_product = flux.operators.wedge_product;
const evolution_mod = flux.integrators.evolution;

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

pub const StateImpl = struct {
    mesh: *const Mesh,
    operators: *operator_context_mod.OperatorContext(Mesh),
    velocity: Velocity,
    vorticity: Vorticity,
    boundary_velocity: []f64,
    velocity_forcing: []f64,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh) !StateImpl {
        var velocity = try Velocity.init(allocator, mesh);
        errdefer velocity.deinit(allocator);

        var vorticity = try Vorticity.init(allocator, mesh);
        errdefer vorticity.deinit(allocator);

        const operators = try operator_context_mod.OperatorContext(Mesh).init(allocator, mesh);
        errdefer operators.deinit();
        _ = try operators.codifferential(2);
        _ = try operators.codifferential(3);
        _ = try operators.exteriorDerivative(flux.forms.Primal, 1);
        _ = try operators.laplacian(1);

        const boundary_velocity = try allocator.alloc(f64, mesh.num_edges());
        errdefer allocator.free(boundary_velocity);
        @memset(boundary_velocity, 0.0);

        const velocity_forcing = try allocator.alloc(f64, mesh.num_edges());
        errdefer allocator.free(velocity_forcing);
        @memset(velocity_forcing, 0.0);

        return .{
            .mesh = mesh,
            .operators = operators,
            .velocity = velocity,
            .vorticity = vorticity,
            .boundary_velocity = boundary_velocity,
            .velocity_forcing = velocity_forcing,
        };
    }

    pub fn deinit(self: *StateImpl, allocator: std.mem.Allocator) void {
        allocator.free(self.velocity_forcing);
        allocator.free(self.boundary_velocity);
        self.vorticity.deinit(allocator);
        self.velocity.deinit(allocator);
        self.operators.deinit();
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

    var state = try StateImpl.init(allocator, &mesh);
    defer state.deinit(allocator);
    try seedReferenceMode(allocator, &state);

    const helicity_initial = try conservedQuantity(allocator, &state);
    const Evolution = evolution_mod.FixedTimeEvolution(*StateImpl, stepImpl);
    var evolution = Evolution.init(&state, config.dt);
    defer evolution.deinit();

    const loop_result = try common.runEvolutionLoop(
        allocator,
        &evolution,
        .{
            .steps = config.steps,
            .dt = config.dt,
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
        Renderer{ .state = &state },
    );

    const helicity_final = try conservedQuantity(allocator, &state);
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

pub fn stepImpl(allocator: std.mem.Allocator, state: *StateImpl, dt: f64) !void {
    _ = dt;
    try recoverVelocityFromVorticity(allocator, state);

    var vorticity = try (try state.operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);

    var advection_density = try wedge_product.wedge(allocator, state.velocity, state.vorticity);
    defer advection_density.deinit(allocator);

    var transport = try (try state.operators.codifferential(3)).apply(allocator, advection_density);
    defer transport.deinit(allocator);
}

pub fn seedReferenceMode(allocator: std.mem.Allocator, state: *StateImpl) !void {
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

    var vorticity = try (try state.operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);
    @memcpy(state.boundary_velocity, state.velocity.values);

    var forcing = try (try state.operators.laplacian(1)).apply(allocator, state.velocity);
    defer forcing.deinit(allocator);
    @memcpy(state.velocity_forcing, forcing.values);
}

pub fn conservedQuantity(allocator: std.mem.Allocator, state: *const StateImpl) !f64 {
    const Helicity = observers.HelicityObserver(StateImpl, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    return observer.evaluate(allocator, state, 0);
}

fn selectVelocity(state: *const StateImpl) *const Velocity {
    return &state.velocity;
}

fn recoverVelocityFromVorticity(allocator: std.mem.Allocator, state: *StateImpl) !void {
    var solve = try poisson.solve_one_form_dirichlet(
        Mesh,
        allocator,
        state.operators,
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
    state: *const StateImpl,

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
