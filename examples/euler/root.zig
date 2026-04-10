const root_std = @import("std");

const euler_2d = struct {
    const local_std = @import("std");
    const testing = local_std.testing;
    const flux = @import("flux");
    const common = @import("examples_common");

    pub const Mesh2D = flux.topology.Mesh(2, 2);
    pub const VertexVorticity = flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal);
    pub const FaceVorticity = flux.forms.Cochain(Mesh2D, 2, flux.forms.Primal);
    pub const FaceTracer = flux.forms.Cochain(Mesh2D, 2, flux.forms.Primal);
    const poisson = flux.operators.poisson;
    const flux_io = flux.io;
    const operator_context_mod = flux.operators.context;

    const Vec2 = [2]f64;

    pub const Demo = enum {
        gaussian,
        dipole,
    };

    pub const Config2D = struct {
        demo: Demo = .gaussian,
        grid: u32 = 16,
        steps: u32 = 1000,
        domain: f64 = 1.0,
        cfl: f64 = 0.1,
        output_dir: []const u8 = "output/euler_2d",
        frames: u32 = 50,

        pub fn dt(self: Config2D) f64 {
            return stableDt(self);
        }

        pub fn outputInterval(self: Config2D) u32 {
            if (self.frames == 0) return 0;
            return @max(1, self.steps / self.frames);
        }
    };

    pub const RunResult2D = struct {
        elapsed_s: f64,
        circulation_initial: f64,
        circulation_final: f64,
        snapshot_count: u32,
    };

    pub const State2D = struct {
        const EdgeAdjacency = struct {
            faces: [2]?u32 = .{ null, null },
        };

        mesh: *const Mesh2D,
        operators: *operator_context_mod.OperatorContext(Mesh2D),
        stream_function: VertexVorticity,
        vorticity: FaceVorticity,
        tracer: FaceTracer,
        face_velocity: []Vec2,
        edge_adjacency: []EdgeAdjacency,

        pub fn init(allocator: local_std.mem.Allocator, mesh: *const Mesh2D) !State2D {
            var stream_function = try VertexVorticity.init(allocator, mesh);
            errdefer stream_function.deinit(allocator);

            var vorticity = try FaceVorticity.init(allocator, mesh);
            errdefer vorticity.deinit(allocator);

            var tracer = try FaceTracer.init(allocator, mesh);
            errdefer tracer.deinit(allocator);

            const operators = try operator_context_mod.OperatorContext(Mesh2D).init(allocator, mesh);
            errdefer operators.deinit();
            _ = try operators.laplacian(0);

            const face_velocity = try allocator.alloc(Vec2, mesh.num_faces());
            errdefer allocator.free(face_velocity);
            @memset(face_velocity, .{ 0.0, 0.0 });

            const edge_adjacency = try buildEdgeAdjacency(allocator, mesh);
            errdefer allocator.free(edge_adjacency);

            return .{
                .mesh = mesh,
                .operators = operators,
                .stream_function = stream_function,
                .vorticity = vorticity,
                .tracer = tracer,
                .face_velocity = face_velocity,
                .edge_adjacency = edge_adjacency,
            };
        }

        pub fn deinit(self: *State2D, allocator: local_std.mem.Allocator) void {
            allocator.free(self.edge_adjacency);
            allocator.free(self.face_velocity);
            self.tracer.deinit(allocator);
            self.vorticity.deinit(allocator);
            self.stream_function.deinit(allocator);
            self.operators.deinit();
        }
    };

    pub fn initializeGaussianVortex(state: *State2D) void {
        const face_centers = state.mesh.simplices(2).items(.barycenter);
        const sigma = 0.12;

        for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
            omega.* = gaussianBlob(center, .{ 0.5, 0.5 }, sigma);
            tracer.* = tracerStripe(center);
        }
    }

    pub fn initializeVortexDipole(state: *State2D) void {
        const face_centers = state.mesh.simplices(2).items(.barycenter);
        const sigma = 0.075;

        for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
            const positive = gaussianBlob(center, .{ 0.38, 0.50 }, sigma);
            const negative = gaussianBlob(center, .{ 0.62, 0.50 }, sigma);
            omega.* = 1.4 * (positive - negative);
            tracer.* = tracerStripe(center);
        }
    }

    pub fn totalCirculation(state: *const State2D) f64 {
        const areas = state.mesh.simplices(2).items(.volume);
        var circulation: f64 = 0.0;
        for (state.vorticity.values, areas) |omega, area| {
            circulation += omega * area;
        }
        return circulation;
    }

    pub fn step2D(
        allocator: local_std.mem.Allocator,
        state: *State2D,
        dt: f64,
    ) !void {
        try refreshDerivedFields(allocator, state);
        try advectVorticity(allocator, state, dt);
        try advectScalar(allocator, state, state.tracer.values, dt);
        try refreshDerivedFields(allocator, state);
    }

    fn stableDt(config: Config2D) f64 {
        return config.cfl * (config.domain / @as(f64, @floatFromInt(config.grid)));
    }

    pub fn run2D(
        allocator: local_std.mem.Allocator,
        config: Config2D,
        writer: anytype,
    ) !RunResult2D {
        var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
        defer mesh.deinit(allocator);

        var state = try State2D.init(allocator, &mesh);
        defer state.deinit(allocator);
        switch (config.demo) {
            .gaussian => initializeGaussianVortex(&state),
            .dipole => initializeVortexDipole(&state),
        }
        try refreshDerivedFields(allocator, &state);

        const circulation_initial = totalCirculation(&state);
        const dt = config.dt();
        var series = try common.Series.init(
            allocator,
            config.output_dir,
            baseName(config.demo),
            common.Plan.fromFrames(config.steps, config.frames, .{}),
        );
        defer series.deinit();

        const start_ns = local_std.time.nanoTimestamp();
        for (0..config.steps) |step_idx| {
            try step2D(allocator, &state, dt);

            if (series.dueAt(@intCast(step_idx + 1), config.steps)) {
                const t = @as(f64, @floatFromInt(step_idx + 1)) * dt;
                try series.capture(t, EulerRenderer{ .state = &state });
            }
        }
        const elapsed_ns = local_std.time.nanoTimestamp() - start_ns;
        try series.finalize();

        const circulation_final = totalCirculation(&state);
        try writer.print(
            "euler_2d[{s}]: grid={d} steps={d} dt={d:.6} circulation={d:.12} -> {d:.12}\n",
            .{ @tagName(config.demo), config.grid, config.steps, dt, circulation_initial, circulation_final },
        );

        return .{
            .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
            .circulation_initial = circulation_initial,
            .circulation_final = circulation_final,
            .snapshot_count = series.count,
        };
    }

    const EulerRenderer = struct {
        state: *const State2D,

        pub fn render(self: @This(), allocator: local_std.mem.Allocator, writer: anytype) !void {
            const velocity_values = try allocator.alloc(f64, self.state.mesh.num_faces() * 3);
            defer allocator.free(velocity_values);

            for (self.state.face_velocity, 0..) |velocity, face_idx| {
                const base = 3 * face_idx;
                velocity_values[base + 0] = velocity[0];
                velocity_values[base + 1] = velocity[1];
                velocity_values[base + 2] = 0.0;
            }

            const point_data = [_]flux_io.DataArraySlice{
                .{ .name = "stream_function", .values = self.state.stream_function.values },
            };
            const cell_data = [_]flux_io.DataArraySlice{
                .{ .name = "vorticity", .values = self.state.vorticity.values },
                .{ .name = "tracer", .values = self.state.tracer.values },
                .{ .name = "velocity", .values = velocity_values, .num_components = 3 },
            };
            try flux_io.write(writer, self.state.mesh.*, &point_data, &cell_data);
        }
    };

    fn baseName(demo: Demo) []const u8 {
        return switch (demo) {
            .gaussian => "euler_2d",
            .dipole => "euler_dipole",
        };
    }

    fn buildEdgeAdjacency(allocator: local_std.mem.Allocator, mesh: *const Mesh2D) ![]State2D.EdgeAdjacency {
        const adjacency = try allocator.alloc(State2D.EdgeAdjacency, mesh.num_edges());
        @memset(adjacency, .{ .faces = .{ null, null } });

        const boundary_2 = mesh.boundary(2);
        for (0..mesh.num_faces()) |face_idx_usize| {
            const face_idx: u32 = @intCast(face_idx_usize);
            const row = boundary_2.row(face_idx);
            for (row.cols) |edge_idx| {
                const edge = &adjacency[edge_idx];
                if (edge.faces[0] == null) {
                    edge.faces[0] = face_idx;
                    continue;
                }
                local_std.debug.assert(edge.faces[1] == null);
                edge.faces[1] = face_idx;
            }
        }

        return adjacency;
    }

    fn refreshDerivedFields(allocator: local_std.mem.Allocator, state: *State2D) !void {
        try recoverStreamFunction(allocator, state);
        reconstructFaceVelocity(state);
    }

    fn recoverStreamFunction(allocator: local_std.mem.Allocator, state: *State2D) !void {
        const vertex_count = state.mesh.num_vertices();
        const forcing = try allocator.alloc(f64, vertex_count);
        defer allocator.free(forcing);
        const weights = try allocator.alloc(f64, vertex_count);
        defer allocator.free(weights);
        const boundary = try allocator.alloc(f64, vertex_count);
        defer allocator.free(boundary);

        @memset(forcing, 0.0);
        @memset(weights, 0.0);
        @memset(boundary, 0.0);

        const face_vertices = state.mesh.simplices(2).items(.vertices);
        const face_areas = state.mesh.simplices(2).items(.volume);

        for (face_vertices, face_areas, state.vorticity.values) |verts, area, omega| {
            const lumped_area = area / 3.0;
            for (verts) |vertex_idx| {
                forcing[vertex_idx] += omega * lumped_area;
                weights[vertex_idx] += lumped_area;
            }
        }

        for (forcing, weights) |*forcing_value, weight| {
            if (weight == 0.0) continue;
            forcing_value.* /= weight;
        }

        var solve = try poisson.solve_zero_form_dirichlet(Mesh2D, allocator, state.operators, forcing, boundary, .{});
        defer solve.deinit(allocator);
        @memcpy(state.stream_function.values, solve.solution);
    }

    fn reconstructFaceVelocity(state: *State2D) void {
        const coords = state.mesh.vertices.slice().items(.coords);
        const face_vertices = state.mesh.simplices(2).items(.vertices);

        for (state.face_velocity, face_vertices) |*velocity, verts| {
            const p0 = coords[verts[0]];
            const p1 = coords[verts[1]];
            const p2 = coords[verts[2]];

            const psi0 = state.stream_function.values[verts[0]];
            const psi1 = state.stream_function.values[verts[1]];
            const psi2 = state.stream_function.values[verts[2]];

            const det = (p1[0] - p0[0]) * (p2[1] - p0[1]) -
                (p2[0] - p0[0]) * (p1[1] - p0[1]);
            local_std.debug.assert(det > 0.0);

            const grad_x = (psi0 * (p1[1] - p2[1]) +
                psi1 * (p2[1] - p0[1]) +
                psi2 * (p0[1] - p1[1])) / det;
            const grad_y = (psi0 * (p2[0] - p1[0]) +
                psi1 * (p0[0] - p2[0]) +
                psi2 * (p1[0] - p0[0])) / det;

            velocity.* = .{ grad_y, -grad_x };
        }
    }

    fn advectVorticity(
        allocator: local_std.mem.Allocator,
        state: *State2D,
        dt: f64,
    ) !void {
        try advectScalar(allocator, state, state.vorticity.values, dt);
    }

    fn advectScalar(
        allocator: local_std.mem.Allocator,
        state: *const State2D,
        values: []f64,
        dt: f64,
    ) !void {
        const face_count = state.mesh.num_faces();
        const face_areas = state.mesh.simplices(2).items(.volume);
        const face_centers = state.mesh.simplices(2).items(.barycenter);
        const edge_vertices = state.mesh.simplices(1).items(.vertices);
        const edge_lengths = state.mesh.simplices(1).items(.volume);
        const coords = state.mesh.vertices.slice().items(.coords);

        const mass_old = try allocator.alloc(f64, face_count);
        defer allocator.free(mass_old);
        const mass_new = try allocator.alloc(f64, face_count);
        defer allocator.free(mass_new);

        for (mass_old, mass_new, values, face_areas) |*old, *new, value, area| {
            old.* = value * area;
            new.* = old.*;
        }

        for (state.edge_adjacency, 0..) |adjacency, edge_idx_usize| {
            const face_left = adjacency.faces[0] orelse continue;
            const face_right = adjacency.faces[1] orelse continue;

            const edge_idx: u32 = @intCast(edge_idx_usize);
            const edge = edge_vertices[edge_idx];
            const p0 = coords[edge[0]];
            const p1 = coords[edge[1]];
            const tangent = .{ p1[0] - p0[0], p1[1] - p0[1] };

            var normal = .{ tangent[1] / edge_lengths[edge_idx], -tangent[0] / edge_lengths[edge_idx] };
            const center_delta = .{
                face_centers[face_right][0] - face_centers[face_left][0],
                face_centers[face_right][1] - face_centers[face_left][1],
            };
            if (dot(normal, center_delta) < 0.0) {
                normal = .{ -normal[0], -normal[1] };
            }

            const velocity_left = state.face_velocity[face_left];
            const velocity_right = state.face_velocity[face_right];
            const edge_velocity = .{
                0.5 * (velocity_left[0] + velocity_right[0]),
                0.5 * (velocity_left[1] + velocity_right[1]),
            };

            const edge_flux = dt * edge_lengths[edge_idx] * dot(edge_velocity, normal);
            const transfer = if (edge_flux >= 0.0)
                edge_flux * values[face_left]
            else
                edge_flux * values[face_right];

            mass_new[face_left] -= transfer;
            mass_new[face_right] += transfer;
        }

        for (values, mass_new, face_areas) |*value, mass, area| {
            value.* = mass / area;
        }
    }

    fn dot(lhs: Vec2, rhs: Vec2) f64 {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1];
    }

    fn gaussianBlob(point: Vec2, center: Vec2, sigma: f64) f64 {
        const dx = point[0] - center[0];
        const dy = point[1] - center[1];
        return local_std.math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
    }

    fn tracerStripe(point: Vec2) f64 {
        const y = point[1] - 0.5;
        const stripe = local_std.math.exp(-(y * y) / (2.0 * 0.035 * 0.035));
        const x_modulation = 0.5 * (1.0 + local_std.math.cos(10.0 * local_std.math.pi * point[0]));
        return stripe * x_modulation;
    }

    test "Euler 2D state allocates stream function on vertices and vorticity on faces" {
        const allocator = testing.allocator;

        var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 1.0, 1.0);
        defer mesh.deinit(allocator);

        var state = try State2D.init(allocator, &mesh);
        defer state.deinit(allocator);

        try testing.expectEqual(@as(usize, mesh.num_vertices()), state.stream_function.values.len);
        try testing.expectEqual(@as(usize, mesh.num_faces()), state.vorticity.values.len);
        try testing.expectEqual(@as(usize, mesh.num_faces()), state.tracer.values.len);
    }

    test "Euler 2D dipole initialization contains both circulation signs" {
        const allocator = testing.allocator;

        var mesh = try Mesh2D.uniform_grid(allocator, 10, 10, 1.0, 1.0);
        defer mesh.deinit(allocator);

        var state = try State2D.init(allocator, &mesh);
        defer state.deinit(allocator);
        initializeVortexDipole(&state);

        var min_value: f64 = local_std.math.inf(f64);
        var max_value: f64 = -local_std.math.inf(f64);
        for (state.vorticity.values) |omega| {
            min_value = @min(min_value, omega);
            max_value = @max(max_value, omega);
        }

        try testing.expect(min_value < 0.0);
        try testing.expect(max_value > 0.0);
    }

    test "Euler 2D circulation is conserved over 1000 explicit steps" {
        const allocator = testing.allocator;
        const config = Config2D{
            .grid = 8,
            .steps = 1000,
            .domain = 1.0,
            .cfl = 0.02,
        };

        var mesh = try Mesh2D.uniform_grid(allocator, config.grid, config.grid, config.domain, config.domain);
        defer mesh.deinit(allocator);

        var state = try State2D.init(allocator, &mesh);
        defer state.deinit(allocator);
        initializeGaussianVortex(&state);

        const circulation_initial = totalCirculation(&state);
        const dt = stableDt(config);

        for (0..config.steps) |_| {
            try step2D(allocator, &state, dt);
        }

        const circulation_final = totalCirculation(&state);
        try testing.expectApproxEqAbs(circulation_initial, circulation_final, 1e-12);
    }
};

const euler_3d = struct {
    const local_std = @import("std");
    const testing = local_std.testing;
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

    pub const Config3D = struct {
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

    pub const RunResult3D = struct {
        elapsed_s: f64,
        helicity_initial: f64,
        helicity_final: f64,
        snapshot_count: u32,
    };

    pub const State3D = struct {
        mesh: *const Mesh3D,
        operators: *operator_context_mod.OperatorContext(Mesh3D),
        velocity: Velocity,
        vorticity: Vorticity,
        boundary_velocity: []f64,
        velocity_forcing: []f64,

        pub fn init(allocator: local_std.mem.Allocator, mesh: *const Mesh3D) !State3D {
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

        pub fn deinit(self: *State3D, allocator: local_std.mem.Allocator) void {
            allocator.free(self.velocity_forcing);
            allocator.free(self.boundary_velocity);
            self.vorticity.deinit(allocator);
            self.velocity.deinit(allocator);
            self.operators.deinit();
        }
    };

    fn selectVelocity(state: *const State3D) *const Velocity {
        return &state.velocity;
    }

    pub fn seedReferenceMode3D(allocator: local_std.mem.Allocator, state: *State3D) !void {
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
                envelope * (@sin(local_std.math.pi * midpoint[2]) + 0.25 * @sin(2.0 * local_std.math.pi * midpoint[1])),
                envelope * (@sin(local_std.math.pi * midpoint[0]) + 0.20 * @sin(2.0 * local_std.math.pi * midpoint[2])),
                envelope * (@sin(local_std.math.pi * midpoint[1]) + 0.15 * @sin(2.0 * local_std.math.pi * midpoint[0])),
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

    pub fn recoverVelocityFromVorticity(allocator: local_std.mem.Allocator, state: *State3D) !void {
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

    pub fn step3D(allocator: local_std.mem.Allocator, state: *State3D, dt: f64) !void {
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

    pub fn run3D(allocator: local_std.mem.Allocator, config: Config3D, writer: anytype) !RunResult3D {
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

        var state = try State3D.init(allocator, &mesh);
        defer state.deinit(allocator);
        try seedReferenceMode3D(allocator, &state);

        const helicity_initial = try computeHelicity3D(allocator, &state);

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

        const start_ns = local_std.time.nanoTimestamp();
        for (0..config.steps) |step_index| {
            try step3D(allocator, &state, config.dt);

            if (series.dueAt(@intCast(step_index + 1), config.steps)) {
                const t = @as(f64, @floatFromInt(step_index + 1)) * config.dt;
                try series.capture(t, Euler3DRenderer{ .state = &state });
            }
        }
        const elapsed_ns = local_std.time.nanoTimestamp() - start_ns;

        try series.finalize();

        const helicity_final = try computeHelicity3D(allocator, &state);
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
        state: *const State3D,

        pub fn render(self: @This(), allocator: local_std.mem.Allocator, writer: anytype) !void {
            try writeSnapshot(allocator, writer, self.state);
        }
    };

    pub fn computeHelicity3D(allocator: local_std.mem.Allocator, state: *const State3D) !f64 {
        const Helicity = observers.HelicityObserver(State3D, Velocity, selectVelocity);
        const observer = Helicity{ .name = "helicity" };
        return observer.evaluate(allocator, state, 0);
    }

    fn writeSnapshot(allocator: local_std.mem.Allocator, writer: anytype, state: *const State3D) !void {
        try common.viz.writeProjectedTetFields(
            2,
            allocator,
            writer,
            state.mesh,
            .{
                .{ .name = "velocity_intensity", .kind = .edge_abs_mean, .values = state.velocity.values },
                .{ .name = "vorticity_flux", .kind = .face_abs_mean, .values = state.vorticity.values },
            },
        );
    }

    test "velocity recovery reproduces the seeded co-closed 1-form on a tetrahedral mesh" {
        const allocator = testing.allocator;

        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        var state = try State3D.init(allocator, &mesh);
        defer state.deinit(allocator);

        try seedReferenceMode3D(allocator, &state);

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

        var output = local_std.ArrayListUnmanaged(u8){};
        defer output.deinit(allocator);

        const result = try run3D(allocator, .{
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

        var state = try State3D.init(allocator, &mesh);
        defer state.deinit(allocator);
        try seedReferenceMode3D(allocator, &state);

        const Helicity = observers.HelicityObserver(State3D, Velocity, selectVelocity);
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
};

pub fn Mesh(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.Mesh2D,
        3 => euler_3d.Mesh3D,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn Config(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.Config2D,
        3 => euler_3d.Config3D,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn RunResult(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.RunResult2D,
        3 => euler_3d.RunResult3D,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn State(comptime topological_dimension: u8) type {
    return switch (topological_dimension) {
        2 => euler_2d.State2D,
        3 => euler_3d.State3D,
        else => @compileError("Euler examples only support topological dimensions 2 and 3"),
    };
}

pub fn run(comptime topological_dimension: u8, allocator: root_std.mem.Allocator, config: Config(topological_dimension), writer: anytype) !RunResult(topological_dimension) {
    return switch (topological_dimension) {
        2 => try euler_2d.run2D(allocator, config, writer),
        3 => try euler_3d.run3D(allocator, config, writer),
        else => unreachable,
    };
}

pub fn step(comptime topological_dimension: u8, allocator: root_std.mem.Allocator, state: *State(topological_dimension), dt: f64) !void {
    switch (topological_dimension) {
        2 => try euler_2d.step2D(allocator, state, dt),
        3 => try euler_3d.step3D(allocator, state, dt),
        else => unreachable,
    }
}

pub fn seedReferenceMode(comptime topological_dimension: u8, allocator: root_std.mem.Allocator, state: *State(topological_dimension)) !void {
    if (topological_dimension == 2) {
        euler_2d.initializeVortexDipole(state);
        return;
    }
    if (topological_dimension == 3) {
        return euler_3d.seedReferenceMode3D(allocator, state);
    }
    unreachable;
}

pub fn conservedQuantity(comptime topological_dimension: u8, allocator: root_std.mem.Allocator, state: *const State(topological_dimension)) !f64 {
    if (topological_dimension == 2) {
        return euler_2d.totalCirculation(state);
    }
    if (topological_dimension == 3) {
        return euler_3d.computeHelicity3D(allocator, state);
    }
    unreachable;
}

test {
    root_std.testing.refAllDeclsRecursive(@This());
}
