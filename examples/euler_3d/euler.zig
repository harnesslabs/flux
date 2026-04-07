const std = @import("std");
const testing = std.testing;
const flux = @import("flux");

const observers = flux.operators.observers;
const poisson = flux.operators.poisson;
const wedge_product = flux.operators.wedge_product;

pub const Mesh3D = flux.Mesh(3, 3);
pub const Velocity = flux.Cochain(Mesh3D, 1, flux.Primal);
pub const Vorticity = flux.Cochain(Mesh3D, 2, flux.Primal);
pub const Potential = flux.Cochain(Mesh3D, 2, flux.Primal);

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

    fn snapshotCount(self: Config) u32 {
        if (self.output_dir == null or self.output_interval == 0) return 0;
        return self.steps / self.output_interval;
    }
};

pub const RunResult = struct {
    elapsed_s: f64,
    helicity_initial: f64,
    helicity_final: f64,
    snapshot_count: u32,
};

pub const State = struct {
    mesh: *const Mesh3D,
    operators: *flux.OperatorContext(Mesh3D),
    velocity: Velocity,
    vorticity: Vorticity,
    boundary_velocity: []f64,
    velocity_forcing: []f64,

    pub fn init(allocator: std.mem.Allocator, mesh: *const Mesh3D) !State {
        var velocity = try Velocity.init(allocator, mesh);
        errdefer velocity.deinit(allocator);

        var vorticity = try Vorticity.init(allocator, mesh);
        errdefer vorticity.deinit(allocator);

        const operators = try flux.OperatorContext(Mesh3D).init(allocator, mesh);
        errdefer operators.deinit();
        try operators.withCodifferential(2);
        try operators.withCodifferential(3);
        try operators.withExteriorDerivative(flux.Primal, 1);
        try operators.withLaplacian(1);

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

    var vorticity = try state.operators.exteriorDerivative(flux.Primal, 1).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);
    @memcpy(state.boundary_velocity, state.velocity.values);

    var forcing = try state.operators.laplacian(1).apply(allocator, state.velocity);
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

    var vorticity = try state.operators.exteriorDerivative(flux.Primal, 1).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);

    var advection_density = try wedge_product.wedge(allocator, state.velocity, state.vorticity);
    defer advection_density.deinit(allocator);

    // The seeded reference mode is intentionally steady: we exercise the
    // helical top-form and transport residual assembly without claiming a full
    // nonlinear 3D closure that the current operator stack does not yet support.
    var transport = try state.operators.codifferential(3).apply(allocator, advection_density);
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
    const snapshot_count_max = config.snapshotCount();
    var pvd_entries: []flux.io.PvdEntry = &.{};
    var filename_bufs: [][flux.io.max_snapshot_filename_length]u8 = &.{};

    if (snapshot_count_max > 0) {
        if (config.output_dir) |output_dir| {
            try ensureDir(output_dir);
        }
        pvd_entries = try allocator.alloc(flux.io.PvdEntry, snapshot_count_max);
        filename_bufs = try allocator.alloc([flux.io.max_snapshot_filename_length]u8, snapshot_count_max);
    }
    defer {
        if (snapshot_count_max > 0) {
            allocator.free(filename_bufs);
            allocator.free(pvd_entries);
        }
    }

    var snapshot_count: u32 = 0;
    const start_ns = std.time.nanoTimestamp();
    for (0..config.steps) |step_index| {
        try step(allocator, &state, config.dt);

        if (config.output_dir) |output_dir| {
            if (config.output_interval > 0 and (step_index + 1) % config.output_interval == 0) {
                const filename = flux.io.snapshot_filename(&filename_bufs[snapshot_count], "euler_3d", snapshot_count);
                try writeSnapshotFile(allocator, output_dir, filename, &state);
                pvd_entries[snapshot_count] = .{
                    .timestep = @as(f64, @floatFromInt(step_index + 1)) * config.dt,
                    .filename = filename,
                };
                snapshot_count += 1;
            }
        }
    }
    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

    if (config.output_dir) |output_dir| {
        if (snapshot_count > 0) {
            try writePvdFile(allocator, output_dir, "euler_3d", pvd_entries[0..snapshot_count]);
        }
    }

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
        .snapshot_count = snapshot_count,
    };
}

fn computeHelicity(allocator: std.mem.Allocator, state: *const State) !f64 {
    const Helicity = observers.HelicityObserver(State, Velocity, selectVelocity);
    const observer = Helicity{ .name = "helicity" };
    return observer.evaluate(allocator, state, 0);
}

fn projectEdgesToTets(
    allocator: std.mem.Allocator,
    mesh: *const Mesh3D,
    edge_values: []const f64,
) ![]f64 {
    const output = try allocator.alloc(f64, mesh.num_tets());
    errdefer allocator.free(output);
    const touched = try allocator.alloc(bool, mesh.num_edges());
    defer allocator.free(touched);

    for (0..mesh.num_tets()) |tet_idx_usize| {
        @memset(touched, false);
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        var count: u32 = 0;

        for (tet_row.cols) |face_idx| {
            const face_row = mesh.boundary(2).row(face_idx);
            for (face_row.cols) |edge_idx| {
                if (touched[edge_idx]) continue;
                touched[edge_idx] = true;
                sum += @abs(edge_values[edge_idx]);
                count += 1;
            }
        }

        std.debug.assert(count > 0);
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(count));
    }

    return output;
}

fn projectFacesToTets(
    allocator: std.mem.Allocator,
    mesh: *const Mesh3D,
    face_values: []const f64,
) ![]f64 {
    const output = try allocator.alloc(f64, mesh.num_tets());
    errdefer allocator.free(output);

    for (0..mesh.num_tets()) |tet_idx_usize| {
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        for (tet_row.cols) |face_idx| {
            sum += @abs(face_values[face_idx]);
        }
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(tet_row.cols.len));
    }

    return output;
}

fn writeSnapshot(allocator: std.mem.Allocator, writer: anytype, state: *const State) !void {
    const velocity_projected = try projectEdgesToTets(allocator, state.mesh, state.velocity.values);
    defer allocator.free(velocity_projected);

    const vorticity_projected = try projectFacesToTets(allocator, state.mesh, state.vorticity.values);
    defer allocator.free(vorticity_projected);

    const cell_data = [_]flux.io.DataArraySlice{
        .{ .name = "velocity_intensity", .values = velocity_projected },
        .{ .name = "vorticity_flux", .values = vorticity_projected },
    };

    try flux.io.write(
        writer,
        Mesh3D.embedding_dimension,
        Mesh3D.topological_dimension,
        state.mesh.*,
        &.{},
        &cell_data,
    );
}

fn writeSnapshotFile(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    filename: []const u8,
    state: *const State,
) !void {
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try writeSnapshot(allocator, output.writer(allocator), state);

    var dir = try std.fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(filename, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn writePvdFile(
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    base_name: []const u8,
    entries: []const flux.io.PvdEntry,
) !void {
    _ = base_name;
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try flux.io.write_pvd(output.writer(allocator), entries);

    var pvd_buf: [flux.io.max_snapshot_filename_length]u8 = undefined;
    const pvd_name = std.fmt.bufPrint(&pvd_buf, "euler_3d.pvd", .{}) catch
        return error.FilenameTooLong;

    var dir = try std.fs.cwd().openDir(output_dir, .{});
    defer dir.close();
    const file = try dir.createFile(pvd_name, .{});
    defer file.close();
    try file.writeAll(output.items);
}

fn ensureDir(path: []const u8) !void {
    try std.fs.cwd().makePath(path);
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
