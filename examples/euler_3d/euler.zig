const std = @import("std");
const testing = std.testing;
const flux = @import("flux");
const common = @import("examples_common");

const observers = flux.operators.observers;
const poisson = flux.operators.poisson;
const wedge_product = flux.operators.wedge_product;
const operator_context_mod = flux.operators.context;

pub const Mesh3D = flux.topology.Mesh(3, 3);
pub const Velocity = flux.forms.Cochain(Mesh3D, 1, flux.forms.Primal);
pub const Vorticity = flux.forms.Cochain(Mesh3D, 2, flux.forms.Primal);
pub const Potential = flux.forms.Cochain(Mesh3D, 2, flux.forms.Primal);

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
};

pub const RunResult = struct {
    elapsed_s: f64,
    helicity_initial: f64,
    helicity_final: f64,
    snapshot_count: u32,
};

pub const State = struct {
    mesh: *const Mesh3D,
    operators: *operator_context_mod.OperatorContext(Mesh3D),
    velocity: Velocity,
    vorticity: Vorticity,
    boundary_velocity: []f64,
    velocity_forcing: []f64,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh3D) !State {
        var velocity = try Velocity.init(allocator, mesh);
        errdefer velocity.deinit(allocator);

        var vorticity = try Vorticity.init(allocator, mesh);
        errdefer vorticity.deinit(allocator);

        const operators = try operator_context_mod.OperatorContext(Mesh3D).init(allocator, mesh);
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

    pub fn deinit(self: *State, allocator: std.mem.Allocator) void {
        allocator.free(self.velocity_forcing);
        allocator.free(self.boundary_velocity);
        self.vorticity.deinit(allocator);
        self.velocity.deinit(allocator);
        self.operators.deinit();
    }
};

fn selectVelocity(state: *const State) *const Velocity {
    return &state.velocity;
}

pub fn seedReferenceMode(allocator: std.mem.Allocator, state: *State) !void {
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

pub fn recoverVelocityFromVorticity(allocator: std.mem.Allocator, state: *State) !void {
    var solve = try poisson.solve_one_form_dirichlet(
        Mesh3D,
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

pub fn step(allocator: std.mem.Allocator, state: *State, dt: f64) !void {
    _ = dt;
    try recoverVelocityFromVorticity(allocator, state);

    var vorticity = try (try state.operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);

    var advection_density = try wedge_product.wedge(allocator, state.velocity, state.vorticity);
    defer advection_density.deinit(allocator);

    // The seeded reference mode is intentionally steady: we exercise the
    // helical top-form and transport residual assembly without claiming a full
    // nonlinear 3D closure that the current operator stack does not yet support.
    var transport = try (try state.operators.codifferential(3)).apply(allocator, advection_density);
    defer transport.deinit(allocator);
}

pub fn run(allocator: std.mem.Allocator, config: Config, writer: anytype) !RunResult {
    var mesh = try Mesh3D.uniform_tetrahedral_grid(
        allocator,
        config.nx,
        config.ny,
        config.nz,
        config.width,
        config.height,
        config.depth,
    );
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    try seedReferenceMode(allocator, &state);

    const helicity_initial = try computeHelicity(allocator, &state);

    const plan: common.Plan = if (config.output_dir != null and config.output_interval > 0)
        common.Plan.fromInterval(config.steps, config.output_interval, .{})
    else
        common.Plan.disabled();

    var series = try common.Series.init(
        allocator,
        config.output_dir orelse "",
        "euler_3d",
        plan,
    );
    defer series.deinit();

    const start_ns = std.time.nanoTimestamp();
    for (0..config.steps) |step_index| {
        try step(allocator, &state, config.dt);

        if (series.dueAt(@intCast(step_index + 1), config.steps)) {
            const t = @as(f64, @floatFromInt(step_index + 1)) * config.dt;
            try series.capture(t, Euler3DRenderer{ .state = &state });
        }
    }
    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

    try series.finalize();

    const helicity_final = try computeHelicity(allocator, &state);
    try writer.print(
        "euler_3d: grid={d}x{d}x{d} steps={d} dt={d:.6} helicity={d:.12} -> {d:.12}\n",
        .{
            config.nx,
            config.ny,
            config.nz,
            config.steps,
            config.dt,
            helicity_initial,
            helicity_final,
        },
    );

    return .{
        .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
        .helicity_initial = helicity_initial,
        .helicity_final = helicity_final,
        .snapshot_count = series.count,
    };
}

const Euler3DRenderer = struct {
    state: *const State,

    pub fn render(self: @This(), allocator: std.mem.Allocator, writer: anytype) !void {
        try writeSnapshot(allocator, writer, self.state);
    }
};

pub fn computeHelicity(allocator: std.mem.Allocator, state: *const State) !f64 {
    const Helicity = observers.HelicityObserver(State, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    return observer.evaluate(allocator, state, 0);
}

fn writeSnapshot(allocator: std.mem.Allocator, writer: anytype, state: *const State) !void {
    const velocity_projected = try common.viz.projectEdgesToTets(allocator, state.mesh, state.velocity.values);
    defer allocator.free(velocity_projected);

    const vorticity_projected = try common.viz.projectFacesToTets(allocator, state.mesh, state.vorticity.values);
    defer allocator.free(vorticity_projected);

    const cell_data = [_]flux.io.DataArraySlice{
        .{ .name = "velocity_intensity", .values = velocity_projected },
        .{ .name = "vorticity_flux", .values = vorticity_projected },
    };

    try flux.io.write(
        writer,
        state.mesh.*,
        &.{},
        &cell_data,
    );
}

test "velocity recovery reproduces the seeded co-closed 1-form on a tetrahedral mesh" {
    const allocator = testing.allocator;

    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);

    try seedReferenceMode(allocator, &state);

    var expected = try Velocity.init(allocator, &mesh);
    defer expected.deinit(allocator);
    @memcpy(expected.values, state.velocity.values);

    try recoverVelocityFromVorticity(allocator, &state);

    for (state.velocity.values, expected.values) |actual, reference| {
        try testing.expectApproxEqAbs(reference, actual, 1e-10);
    }
}

test "helicity is conserved over 1000 steps for the seeded reference mode" {
    const allocator = testing.allocator;

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const result = try run(allocator, .{
        .steps = 1000,
        .nx = 2,
        .ny = 2,
        .nz = 2,
        .dt = 0.01,
        .output_dir = null,
        .output_interval = 0,
    }, output.writer(allocator));

    try testing.expectEqual(@as(u32, 0), result.snapshot_count);
    try testing.expectApproxEqAbs(result.helicity_initial, result.helicity_final, 1e-12);
}

test "seeded reference mode has nonzero helicity and matches the observer" {
    const allocator = testing.allocator;

    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try State.init(allocator, &mesh);
    defer state.deinit(allocator);
    try seedReferenceMode(allocator, &state);

    const Helicity = observers.HelicityObserver(State, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    const observed = try observer.evaluate(allocator, &state, 0);

    var density = try wedge_product.wedge(allocator, state.velocity, state.vorticity);
    defer density.deinit(allocator);

    var manual: f64 = 0.0;
    for (density.values) |value| {
        manual += value;
    }

    try testing.expect(@abs(manual) > 1e-9);
    try testing.expectApproxEqAbs(manual, observed, 1e-12);
}
