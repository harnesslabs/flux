const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const cochain = flux.forms;
const topology = flux.topology;
const dec_context_mod = flux.operators.dec.context;
const feec_context_mod = flux.operators.feec.context;

pub const Mesh2D = topology.Mesh(2, 2);
pub const Mesh3D = topology.Mesh(3, 3);

pub const BoundaryCondition = enum {
    pec,
};

fn MaxwellCore(comptime MeshType: type, comptime Hooks: type) type {
    return struct {
        pub const State = struct {
            const Self = @This();

            pub const OneForm = cochain.Cochain(MeshType, 1, cochain.Primal);
            pub const TwoForm = cochain.Cochain(MeshType, 2, cochain.Primal);

            E: OneForm,
            B: TwoForm,
            J: OneForm,
            mesh: *const MeshType,
            dec_operators: *dec_context_mod.OperatorContext(MeshType),
            feec_operators: *feec_context_mod.OperatorContext(MeshType),
            timestep: u64,

            pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
                var electric = try OneForm.init(allocator, mesh);
                errdefer electric.deinit(allocator);

                var magnetic = try TwoForm.init(allocator, mesh);
                errdefer magnetic.deinit(allocator);

                var current = try OneForm.init(allocator, mesh);
                errdefer current.deinit(allocator);

                const dec_operators = try dec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
                errdefer dec_operators.deinit();
                const feec_operators = try feec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
                errdefer feec_operators.deinit();
                try Hooks.primeOperators(dec_operators, feec_operators);

                return .{
                    .E = electric,
                    .B = magnetic,
                    .J = current,
                    .mesh = mesh,
                    .dec_operators = dec_operators,
                    .feec_operators = feec_operators,
                    .timestep = 0,
                };
            }

            pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
                self.feec_operators.deinit();
                self.dec_operators.deinit();
                self.J.deinit(allocator);
                self.B.deinit(allocator);
                self.E.deinit(allocator);
            }
        };

        pub fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var derivative = try (try state.dec_operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, state.E);
            defer derivative.deinit(allocator);
            derivative.scale(dt);
            state.B.sub(derivative);
        }

        pub fn applyPecBoundary(state: anytype) void {
            for (state.mesh.boundary_edges) |edge_idx| {
                state.E.values[edge_idx] = 0.0;
            }
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            try Hooks.ampereStep(allocator, state, dt);
        }

        pub fn leapfrogStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            try faradayStep(allocator, state, dt);
            try ampereStep(allocator, state, dt);
            if (Hooks.apply_boundary_in_leapfrog) applyPecBoundary(state);
            state.timestep += 1;
        }
    };
}

fn hooks2d(comptime MeshType: type) type {
    return struct {
        pub const apply_boundary_in_leapfrog = false;

        pub fn primeOperators(
            dec_operators: *dec_context_mod.OperatorContext(MeshType),
            feec_operators: *feec_context_mod.OperatorContext(MeshType),
        ) !void {
            _ = try dec_operators.exteriorDerivative(cochain.Primal, 1);
            _ = try dec_operators.exteriorDerivative(cochain.Dual, 0);
            _ = try feec_operators.hodgeStar(1);
            _ = try feec_operators.hodgeStar(2);
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var star_B = try (try state.feec_operators.hodgeStar(2)).apply(allocator, state.B);
            defer star_B.deinit(allocator);

            var d_star_B = try (try state.dec_operators.exteriorDerivative(cochain.Dual, 0)).apply(allocator, star_B);
            defer d_star_B.deinit(allocator);

            const edge_volumes = state.mesh.simplices(1).items(.volume);
            const dual_edge_volumes = state.mesh.dual_edge_volumes;
            for (state.E.values, d_star_B.values, edge_volumes, dual_edge_volumes, state.J.values) |*e, dsb, volume, dual_volume, j| {
                std.debug.assert(dual_volume > 0.0);
                e.* += dt * ((volume / dual_volume) * dsb - j);
            }
        }
    };
}

fn hooks3d(comptime MeshType: type) type {
    return struct {
        pub const apply_boundary_in_leapfrog = true;

        pub fn primeOperators(
            dec_operators: *dec_context_mod.OperatorContext(MeshType),
            feec_operators: *feec_context_mod.OperatorContext(MeshType),
        ) !void {
            _ = try dec_operators.exteriorDerivative(cochain.Primal, 1);
            _ = try dec_operators.exteriorDerivative(cochain.Primal, 2);
            _ = try dec_operators.exteriorDerivative(cochain.Dual, 1);
            _ = try feec_operators.hodgeStar(2);
            _ = try feec_operators.hodgeStarInverse(1);
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var star_b = try (try state.feec_operators.hodgeStar(2)).apply(allocator, state.B);
            defer star_b.deinit(allocator);

            var derivative = try (try state.dec_operators.exteriorDerivative(cochain.Dual, 1)).apply(allocator, star_b);
            defer derivative.deinit(allocator);

            var curl_b = try (try state.feec_operators.hodgeStarInverse(1)).apply(allocator, derivative);
            defer curl_b.deinit(allocator);

            for (state.E.values, curl_b.values, state.J.values) |*electric, curl_value, current| {
                electric.* += dt * (curl_value - current);
            }
        }
    };
}

pub fn StateForMesh2D(comptime MeshType: type) type {
    return MaxwellCore(MeshType, hooks2d(MeshType)).State;
}

pub fn StateForMesh3D(comptime MeshType: type) type {
    return MaxwellCore(MeshType, hooks3d(MeshType)).State;
}

pub const MaxwellState2D = StateForMesh2D(Mesh2D);
pub const MaxwellState3D = StateForMesh3D(Mesh3D);

pub fn faraday_step(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).faradayStep(allocator, state, dt);
}

pub fn ampere_step(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).ampereStep(allocator, state, dt);
}

pub fn ampere_step_3d(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks3d(@TypeOf(state.mesh.*))).ampereStep(allocator, state, dt);
}

pub fn apply_pec_boundary(state: anytype) void {
    MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).applyPecBoundary(state);
}

pub fn leapfrog_step(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).leapfrogStep(allocator, state, dt);
}

pub fn leapfrog_step_3d(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks3d(@TypeOf(state.mesh.*))).leapfrogStep(allocator, state, dt);
}

pub fn electromagnetic_energy(allocator: std.mem.Allocator, state: anytype) !f64 {
    var star_E = try (try state.feec_operators.hodgeStar(1)).apply(allocator, state.E);
    defer star_E.deinit(allocator);
    var star_B = try (try state.feec_operators.hodgeStar(2)).apply(allocator, state.B);
    defer star_B.deinit(allocator);

    var e_energy: f64 = 0.0;
    for (state.E.values, star_E.values) |e, se| e_energy += e * se;

    var b_energy: f64 = 0.0;
    for (state.B.values, star_B.values) |b, sb| b_energy += b * sb;

    return 0.5 * (e_energy + b_energy);
}

pub fn PointDipole(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        edge_index: u32,
        edge_length: f64,
        frequency: f64,
        amplitude: f64,

        pub fn init(mesh: *const MeshType, frequency: f64, amplitude: f64, position: [MeshType.embedding_dimension]f64) Self {
            const simplex_1 = mesh.simplices(1);
            const edge_verts = simplex_1.items(.vertices);
            const lengths = simplex_1.items(.volume);
            const coords = mesh.vertices.slice().items(.coords);

            var best_edge: u32 = 0;
            var best_distance_squared: f64 = std.math.inf(f64);
            for (0..mesh.num_edges()) |e| {
                const v0 = coords[edge_verts[e][0]];
                const v1 = coords[edge_verts[e][1]];
                var midpoint: [MeshType.embedding_dimension]f64 = undefined;
                inline for (0..MeshType.embedding_dimension) |d| midpoint[d] = 0.5 * (v0[d] + v1[d]);
                var dist_sq: f64 = 0.0;
                inline for (0..MeshType.embedding_dimension) |d| {
                    const diff = midpoint[d] - position[d];
                    dist_sq += diff * diff;
                }
                if (dist_sq < best_distance_squared) {
                    best_distance_squared = dist_sq;
                    best_edge = @intCast(e);
                }
            }

            return .{
                .edge_index = best_edge,
                .edge_length = lengths[best_edge],
                .frequency = frequency,
                .amplitude = amplitude,
            };
        }

        pub fn apply(self: Self, J: anytype, t: f64) void {
            @memset(J.values, 0.0);
            J.values[self.edge_index] = (self.amplitude / self.edge_length) * @sin(2.0 * std.math.pi * self.frequency * t);
        }
    };
}

pub fn DipoleOptions(comptime dim: u8) type {
    return struct {
        center: [dim]f64,
        frequency: f64,
        amplitude: f64,
        boundary: BoundaryCondition = .pec,
    };
}

pub fn CavityOptions(comptime dim: u8) type {
    return switch (dim) {
        2 => struct {
            domain_length: f64,
            time_step: f64,
            boundary: BoundaryCondition = .pec,
        },
        3 => struct {
            width: f64,
            height: f64,
            time_step: f64,
            boundary: BoundaryCondition = .pec,
        },
        else => @compileError("Maxwell examples only support topological dimensions 2 and 3"),
    };
}

fn Maxwell2D(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const Field = enum { electric, magnetic, current };
        pub const Measurement = enum { energy };

        state: StateForMesh2D(MeshType),
        source: ?PointDipole(MeshType) = null,
        boundary: BoundaryCondition = .pec,

        pub const Leapfrog = struct {
            pub fn first(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                const time = @as(f64, @floatFromInt(system.state.timestep)) * dt;
                if (system.source) |source| {
                    source.apply(&system.state.J, time);
                } else {
                    @memset(system.state.J.values, 0.0);
                }
                try faraday_step(allocator, &system.state, dt);
            }

            pub fn second(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try ampere_step(allocator, &system.state, dt);
                system.state.timestep += 1;
            }

            pub fn applyBoundary(system: *Self) void {
                switch (system.boundary) {
                    .pec => apply_pec_boundary(&system.state),
                }
            }
        };

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = try StateForMesh2D(MeshType).init(allocator, mesh),
            };
        }

        pub fn dipole(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: DipoleOptions(2),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            system.source = PointDipole(MeshType).init(mesh, options.frequency, options.amplitude, options.center);
            system.boundary = options.boundary;
            return system;
        }

        pub fn cavity(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: CavityOptions(2),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            @import("reference.zig").project_te10_b(mesh, system.state.B.values, -options.time_step / 2.0, options.domain_length);
            system.boundary = options.boundary;
            system.source = null;
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return switch (which) {
                .energy => try electromagnetic_energy(allocator, &self.state),
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

        state: StateForMesh3D(MeshType),
        boundary: BoundaryCondition = .pec,

        pub const Leapfrog = struct {
            pub fn first(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try faraday_step(allocator, &system.state, dt);
            }

            pub fn second(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try ampere_step_3d(allocator, &system.state, dt);
                system.state.timestep += 1;
            }

            pub fn applyBoundary(system: *Self) void {
                switch (system.boundary) {
                    .pec => apply_pec_boundary(&system.state),
                }
            }
        };

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = try StateForMesh3D(MeshType).init(allocator, mesh),
            };
        }

        pub fn cavity(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: CavityOptions(3),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            try @import("reference.zig").seedTm110Mode(allocator, &system.state, options.time_step, options.width, options.height);
            system.boundary = options.boundary;
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return switch (which) {
                .energy => try electromagnetic_energy(allocator, &self.state),
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
