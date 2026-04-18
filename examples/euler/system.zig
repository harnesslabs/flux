const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const bridges = flux.operators.bridges;
const dec_context_mod = flux.operators.dec.context;
const feec_context_mod = flux.operators.feec.context;
const observers = flux.operators.observers;
const poisson = flux.operators.poisson;
const wedge_product = flux.operators.wedge_product;

pub fn Euler(comptime dim: u8, comptime MeshType: type) type {
    comptime {
        if (dim != 2 and dim != 3) @compileError("euler examples only support dimensions 2 and 3");
    }

    return struct {
        const Self = @This();

        pub const State = if (dim == 2) State2D(MeshType) else State3D(MeshType);
        pub const Field = if (dim == 2)
            enum { stream_function, vorticity, tracer, velocity }
        else
            enum { velocity, vorticity };
        pub const Measurement = if (dim == 2)
            enum { circulation }
        else
            enum { helicity };

        state: State,

        pub const Explicit = struct {
            pub fn advance(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try system.step(allocator, dt);
            }
        };

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = if (dim == 2)
                    try State2D(MeshType).init(allocator, mesh)
                else
                    try State3D(MeshType).init(allocator, mesh),
            };
        }

        pub fn gaussian(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            comptime if (dim != 2) @compileError("gaussian scenario only applies to 2D Euler");

            var system = try Self.init(allocator, mesh);
            initializeGaussianVortex(&system.state);
            try refreshDerivedFields2D(allocator, &system.state);
            return system;
        }

        pub fn dipole(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            comptime if (dim != 2) @compileError("dipole scenario only applies to 2D Euler");

            var system = try Self.init(allocator, mesh);
            initializeVortexDipole(&system.state);
            try refreshDerivedFields2D(allocator, &system.state);
            return system;
        }

        pub fn reference(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            comptime if (dim != 3) @compileError("reference scenario only applies to 3D Euler");

            var system = try Self.init(allocator, mesh);
            try seedReferenceMode3D(allocator, &system.state);
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return if (dim == 2)
                switch (which) {
                    .circulation => totalCirculation(&self.state),
                }
            else switch (which) {
                .helicity => try totalHelicity(allocator, &self.state),
            };
        }

        pub fn writeFields(self: *const Self, allocator: std.mem.Allocator, writer: anytype, fields: []const Field) !void {
            if (dim == 2) {
                try writeFields2D(&self.state, allocator, writer, fields);
                return;
            }
            try writeFields3D(&self.state, allocator, writer, fields);
        }

        fn step(self: *Self, allocator: std.mem.Allocator, dt: f64) !void {
            if (dim == 2) {
                try step2D(allocator, &self.state, dt);
                return;
            }
            try step3D(allocator, &self.state, dt);
        }
    };
}

fn State2D(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const StreamFunction = flux.forms.Cochain(MeshType, 0, flux.forms.Primal);
        pub const Vorticity = flux.forms.Cochain(MeshType, 2, flux.forms.Primal);
        pub const Tracer = flux.forms.Cochain(MeshType, 2, flux.forms.Primal);
        pub const Vec2 = [2]f64;
        const EdgeAdjacency = struct {
            faces: [2]?u32 = .{ null, null },
        };

        mesh: *const MeshType,
        feec_operators: *feec_context_mod.OperatorContext(MeshType),
        stream_function: StreamFunction,
        vorticity: Vorticity,
        tracer: Tracer,
        face_velocity: []Vec2,
        edge_adjacency: []EdgeAdjacency,

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var stream_function = try StreamFunction.init(allocator, mesh);
            errdefer stream_function.deinit(allocator);

            var vorticity = try Vorticity.init(allocator, mesh);
            errdefer vorticity.deinit(allocator);

            var tracer = try Tracer.init(allocator, mesh);
            errdefer tracer.deinit(allocator);

            const feec_operators = try feec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
            errdefer feec_operators.deinit();
            _ = try feec_operators.laplacian(0);

            const face_velocity = try allocator.alloc(Vec2, mesh.num_faces());
            errdefer allocator.free(face_velocity);
            @memset(face_velocity, .{ 0.0, 0.0 });

            const edge_adjacency = try buildEdgeAdjacency(allocator, mesh);
            errdefer allocator.free(edge_adjacency);

            return .{
                .mesh = mesh,
                .feec_operators = feec_operators,
                .stream_function = stream_function,
                .vorticity = vorticity,
                .tracer = tracer,
                .face_velocity = face_velocity,
                .edge_adjacency = edge_adjacency,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.edge_adjacency);
            allocator.free(self.face_velocity);
            self.tracer.deinit(allocator);
            self.vorticity.deinit(allocator);
            self.stream_function.deinit(allocator);
            self.feec_operators.deinit();
        }
    };
}

fn State3D(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const Velocity = flux.forms.Cochain(MeshType, 1, flux.forms.Primal);
        pub const Vorticity = flux.forms.Cochain(MeshType, 2, flux.forms.Primal);

        mesh: *const MeshType,
        dec_operators: *dec_context_mod.OperatorContext(MeshType),
        feec_operators: *feec_context_mod.OperatorContext(MeshType),
        velocity: Velocity,
        vorticity: Vorticity,
        boundary_velocity: []f64,
        velocity_forcing: []f64,

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var velocity = try Velocity.init(allocator, mesh);
            errdefer velocity.deinit(allocator);

            var vorticity = try Vorticity.init(allocator, mesh);
            errdefer vorticity.deinit(allocator);

            const dec_operators = try dec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
            errdefer dec_operators.deinit();
            const feec_operators = try feec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
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

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.velocity_forcing);
            allocator.free(self.boundary_velocity);
            self.vorticity.deinit(allocator);
            self.velocity.deinit(allocator);
            self.feec_operators.deinit();
            self.dec_operators.deinit();
        }
    };
}

fn writeFields2D(state: anytype, allocator: std.mem.Allocator, writer: anytype, fields: anytype) !void {
    var point_data: [1]flux.io.DataArraySlice = undefined;
    var cell_data: [3]flux.io.DataArraySlice = undefined;
    var point_count: usize = 0;
    var cell_count: usize = 0;

    var velocity_values: ?[]f64 = null;
    defer if (velocity_values) |values| allocator.free(values);

    for (fields) |field| {
        switch (field) {
            .stream_function => {
                point_data[point_count] = .{
                    .name = "stream_function",
                    .values = state.stream_function.values,
                };
                point_count += 1;
            },
            .vorticity => {
                cell_data[cell_count] = .{
                    .name = "vorticity",
                    .values = state.vorticity.values,
                };
                cell_count += 1;
            },
            .tracer => {
                cell_data[cell_count] = .{
                    .name = "tracer",
                    .values = state.tracer.values,
                };
                cell_count += 1;
            },
            .velocity => {
                if (velocity_values == null) {
                    velocity_values = try allocator.alloc(f64, state.mesh.num_faces() * 3);
                    for (state.face_velocity, 0..) |velocity, face_idx| {
                        const base = 3 * face_idx;
                        velocity_values.?[base + 0] = velocity[0];
                        velocity_values.?[base + 1] = velocity[1];
                        velocity_values.?[base + 2] = 0.0;
                    }
                }
                cell_data[cell_count] = .{
                    .name = "velocity",
                    .values = velocity_values.?,
                    .num_components = 3,
                };
                cell_count += 1;
            },
        }
    }

    try flux.io.write(writer, state.mesh.*, point_data[0..point_count], cell_data[0..cell_count]);
}

fn writeFields3D(state: anytype, allocator: std.mem.Allocator, writer: anytype, fields: anytype) !void {
    var projected_fields: [2]common.viz.TetProjectionField = undefined;
    var field_count: usize = 0;
    for (fields) |field| {
        projected_fields[field_count] = switch (field) {
            .velocity => .{ .name = "velocity_intensity", .kind = .edge_abs_mean, .values = state.velocity.values },
            .vorticity => .{ .name = "vorticity_flux", .kind = .face_abs_mean, .values = state.vorticity.values },
        };
        field_count += 1;
    }

    switch (field_count) {
        0 => {
            try flux.io.write(writer, state.mesh.*, &.{}, &.{});
        },
        1 => try common.viz.writeProjectedTetFields(1, allocator, writer, state.mesh, .{
            projected_fields[0],
        }),
        2 => try common.viz.writeProjectedTetFields(2, allocator, writer, state.mesh, .{
            projected_fields[0],
            projected_fields[1],
        }),
        else => unreachable,
    }
}

fn buildEdgeAdjacency(allocator: std.mem.Allocator, mesh: anytype) ![]State2D(@TypeOf(mesh.*)).EdgeAdjacency {
    const adjacency = try allocator.alloc(State2D(@TypeOf(mesh.*)).EdgeAdjacency, mesh.num_edges());
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
            std.debug.assert(edge.faces[1] == null);
            edge.faces[1] = face_idx;
        }
    }
    return adjacency;
}

fn refreshDerivedFields2D(allocator: std.mem.Allocator, state: anytype) !void {
    try recoverStreamFunction2D(allocator, state);
    reconstructFaceVelocity2D(state);
}

fn recoverStreamFunction2D(allocator: std.mem.Allocator, state: anytype) !void {
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

    var solve = try poisson.solve_zero_form_dirichlet(
        @TypeOf(state.mesh.*),
        allocator,
        state.feec_operators,
        forcing,
        boundary,
        .{},
    );
    defer solve.deinit(allocator);
    @memcpy(state.stream_function.values, solve.solution);
}

fn reconstructFaceVelocity2D(state: anytype) void {
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
        std.debug.assert(det > 0.0);

        const grad_x = (psi0 * (p1[1] - p2[1]) +
            psi1 * (p2[1] - p0[1]) +
            psi2 * (p0[1] - p1[1])) / det;
        const grad_y = (psi0 * (p2[0] - p1[0]) +
            psi1 * (p0[0] - p2[0]) +
            psi2 * (p1[0] - p0[0])) / det;

        velocity.* = .{ grad_y, -grad_x };
    }
}

fn step2D(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try refreshDerivedFields2D(allocator, state);
    try advectScalar2D(allocator, state, state.vorticity.values, dt);
    try advectScalar2D(allocator, state, state.tracer.values, dt);
    try refreshDerivedFields2D(allocator, state);
}

fn advectScalar2D(allocator: std.mem.Allocator, state: anytype, values: []f64, dt: f64) !void {
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
        if (dot2(normal, center_delta) < 0.0) {
            normal = .{ -normal[0], -normal[1] };
        }

        const velocity_left = state.face_velocity[face_left];
        const velocity_right = state.face_velocity[face_right];
        const edge_velocity = .{
            0.5 * (velocity_left[0] + velocity_right[0]),
            0.5 * (velocity_left[1] + velocity_right[1]),
        };

        const edge_flux = dt * edge_lengths[edge_idx] * dot2(edge_velocity, normal);
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

fn initializeGaussianVortex(state: anytype) void {
    const face_centers = state.mesh.simplices(2).items(.barycenter);
    const sigma = 0.12;

    for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
        omega.* = gaussianBlob(center, .{ 0.5, 0.5 }, sigma);
        tracer.* = tracerStripe(center);
    }
}

fn initializeVortexDipole(state: anytype) void {
    const face_centers = state.mesh.simplices(2).items(.barycenter);
    const sigma = 0.075;

    for (state.vorticity.values, state.tracer.values, face_centers) |*omega, *tracer, center| {
        const positive = gaussianBlob(center, .{ 0.38, 0.50 }, sigma);
        const negative = gaussianBlob(center, .{ 0.62, 0.50 }, sigma);
        omega.* = 1.4 * (positive - negative);
        tracer.* = tracerStripe(center);
    }
}

fn totalCirculation(state: anytype) f64 {
    const areas = state.mesh.simplices(2).items(.volume);
    var circulation: f64 = 0.0;
    for (state.vorticity.values, areas) |omega, area| {
        circulation += omega * area;
    }
    return circulation;
}

fn step3D(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    _ = dt;
    try recoverVelocityFromVorticity3D(allocator, state);

    var vorticity = try (try state.dec_operators.exteriorDerivative(flux.forms.Primal, 1)).apply(allocator, state.velocity);
    defer vorticity.deinit(allocator);
    @memcpy(state.vorticity.values, vorticity.values);

    const interpolate_velocity = bridges.WhitneyInterpolation(@TypeOf(state.mesh.*), 1).init(state.mesh);
    const interpolate_vorticity = bridges.WhitneyInterpolation(@TypeOf(state.mesh.*), 2).init(state.mesh);
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

fn seedReferenceMode3D(allocator: std.mem.Allocator, state: anytype) !void {
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

fn recoverVelocityFromVorticity3D(allocator: std.mem.Allocator, state: anytype) !void {
    var solve = try poisson.solve_one_form_dirichlet(
        @TypeOf(state.mesh.*),
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

fn totalHelicity(allocator: std.mem.Allocator, state: anytype) !f64 {
    const Velocity = @TypeOf(state.velocity);
    const Helicity = observers.HelicityObserver(@TypeOf(state.*), Velocity, selectVelocity3D(@TypeOf(state.*), Velocity));
    const observer = Helicity{ .name = "helicity" };
    return observer.evaluate(allocator, state, 0);
}

fn selectVelocity3D(comptime StateType: type, comptime VelocityType: type) fn (*const StateType) *const VelocityType {
    return struct {
        fn apply(state: *const StateType) *const VelocityType {
            return &state.velocity;
        }
    }.apply;
}

fn dot2(lhs: [2]f64, rhs: [2]f64) f64 {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1];
}

fn gaussianBlob(point: [2]f64, center: [2]f64, sigma: f64) f64 {
    const dx = point[0] - center[0];
    const dy = point[1] - center[1];
    return std.math.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
}

fn tracerStripe(point: [2]f64) f64 {
    const y = point[1] - 0.5;
    const stripe = std.math.exp(-(y * y) / (2.0 * 0.035 * 0.035));
    const x_modulation = 0.5 * (1.0 + std.math.cos(10.0 * std.math.pi * point[0]));
    return stripe * x_modulation;
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
