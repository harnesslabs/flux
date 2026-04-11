const std = @import("std");
const flux = @import("flux");

const cochain = flux.forms;
const topology = flux.topology;
const operator_context_mod = flux.operators.context;

pub const Mesh2D = topology.Mesh(2, 2);
pub const Mesh3D = topology.Mesh(3, 3);

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
            operators: *operator_context_mod.OperatorContext(MeshType),
            timestep: u64,

            pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
                var electric = try OneForm.init(allocator, mesh);
                errdefer electric.deinit(allocator);

                var magnetic = try TwoForm.init(allocator, mesh);
                errdefer magnetic.deinit(allocator);

                var current = try OneForm.init(allocator, mesh);
                errdefer current.deinit(allocator);

                const operators = try operator_context_mod.OperatorContext(MeshType).init(allocator, mesh);
                errdefer operators.deinit();
                try Hooks.primeOperators(operators);

                return .{
                    .E = electric,
                    .B = magnetic,
                    .J = current,
                    .mesh = mesh,
                    .operators = operators,
                    .timestep = 0,
                };
            }

            pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
                self.operators.deinit();
                self.J.deinit(allocator);
                self.B.deinit(allocator);
                self.E.deinit(allocator);
            }
        };

        pub fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var derivative = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, state.E);
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

        pub fn primeOperators(operators: *operator_context_mod.OperatorContext(MeshType)) !void {
            _ = try operators.exteriorDerivative(cochain.Primal, 1);
            _ = try operators.exteriorDerivative(cochain.Dual, 0);
            _ = try operators.hodgeStar(1);
            _ = try operators.hodgeStar(2);
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var star_B = try (try state.operators.hodgeStar(2)).apply(allocator, state.B);
            defer star_B.deinit(allocator);

            var d_star_B = try (try state.operators.exteriorDerivative(cochain.Dual, 0)).apply(allocator, star_B);
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

pub fn StateForMesh2D(comptime MeshType: type) type {
    return MaxwellCore(MeshType, hooks2d(MeshType)).State;
}

pub fn StateForMesh3D(comptime MeshType: type) type {
    return MaxwellCore(MeshType, hooks3d(MeshType)).State;
}

pub const MaxwellState2D = StateForMesh2D(Mesh2D);
pub const MaxwellState3D = StateForMesh3D(Mesh3D);

pub const DrivenLeapfrog2D = struct {
    state: *MaxwellState2D,
    source: ?PointDipole(Mesh2D),
    dt: f64,
    next_step_index: u32 = 0,

    pub fn step(self: *@This(), allocator: std.mem.Allocator) !void {
        const time = @as(f64, @floatFromInt(self.next_step_index)) * self.dt;
        if (self.source) |dipole| dipole.apply(&self.state.J, time);
        try leapfrog_step(allocator, self.state, self.dt);
        apply_pec_boundary(self.state);
        self.next_step_index += 1;
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        _ = self;
        _ = allocator;
    }
};

pub fn faraday_step(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).faradayStep(allocator, state, dt);
}

pub fn ampere_step(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
    try MaxwellCore(@TypeOf(state.mesh.*), hooks2d(@TypeOf(state.mesh.*))).ampereStep(allocator, state, dt);
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

pub fn leapfrog_step_3d_typed(allocator: std.mem.Allocator, state: *MaxwellState3D, dt: f64) !void {
    try leapfrog_step_3d(allocator, state, dt);
}

pub fn electromagnetic_energy(allocator: std.mem.Allocator, state: anytype) !f64 {
    var star_E = try (try state.operators.hodgeStar(1)).apply(allocator, state.E);
    defer star_E.deinit(allocator);
    var star_B = try (try state.operators.hodgeStar(2)).apply(allocator, state.B);
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
