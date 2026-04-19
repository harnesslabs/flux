const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");

const cochain = flux.forms;
const bridges = flux.operators.bridges;
const operator_context_mod = flux.operators.context;

pub const BoundaryCondition = enum {
    pec,
};

pub fn Maxwell(comptime dim: u8, comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const State = struct {
            const StateSelf = @This();

            pub const OneForm = cochain.Cochain(MeshType, 1, cochain.Primal);
            pub const TwoForm = cochain.Cochain(MeshType, 2, cochain.Primal);

            E: OneForm,
            B: TwoForm,
            J: OneForm,
            mesh: *const MeshType,
            operators: *operator_context_mod.OperatorContext(MeshType),
            timestep: u64,

            pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !StateSelf {
                var electric = try OneForm.init(allocator, mesh);
                errdefer electric.deinit(allocator);

                var magnetic = try TwoForm.init(allocator, mesh);
                errdefer magnetic.deinit(allocator);

                var current = try OneForm.init(allocator, mesh);
                errdefer current.deinit(allocator);

                const operators = try operator_context_mod.OperatorContext(MeshType).init(allocator, mesh);
                errdefer operators.deinit();
                try primeOperators(MeshType, operators);

                return .{
                    .E = electric,
                    .B = magnetic,
                    .J = current,
                    .mesh = mesh,
                    .operators = operators,
                    .timestep = 0,
                };
            }

            pub fn deinit(self: *StateSelf, allocator: std.mem.Allocator) void {
                self.operators.deinit();
                self.J.deinit(allocator);
                self.B.deinit(allocator);
                self.E.deinit(allocator);
            }
        };

        pub const Field = enum { electric, magnetic, current };
        pub const Measurement = enum { energy };

        state: State,
        source: if (dim == 2) ?PointDipole(MeshType) else void = if (dim == 2) null else {},
        boundary: BoundaryCondition = .pec,

        pub const Leapfrog = struct {
            pub fn first(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                if (dim == 2) {
                    const time = @as(f64, @floatFromInt(system.state.timestep)) * dt;
                    if (system.source) |source| {
                        source.apply(&system.state.J, time);
                    } else {
                        @memset(system.state.J.values, 0.0);
                    }
                } else {
                    @memset(system.state.J.values, 0.0);
                }
                try Self.faradayStep(allocator, &system.state, dt);
            }

            pub fn second(allocator: std.mem.Allocator, system: *Self, dt: f64) !void {
                try Self.ampereStep(allocator, &system.state, dt);
                system.state.timestep += 1;
            }

            pub fn applyBoundary(system: *Self) void {
                switch (system.boundary) {
                    .pec => Self.applyPecBoundary(&system.state),
                }
            }
        };

        fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var derivative = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, state.E);
            defer derivative.deinit(allocator);
            derivative.scale(dt);
            state.B.sub(derivative);
        }

        fn applyPecBoundary(state: anytype) void {
            for (state.mesh.boundary_edges) |edge_idx| {
                state.E.values[edge_idx] = 0.0;
            }
        }

        fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var star_b = try (try state.operators.hodgeStar(2)).apply(allocator, state.B);
            defer star_b.deinit(allocator);

            switch (MeshType.topological_dimension) {
                2 => {
                    var d_star_b = try (try state.operators.exteriorDerivative(cochain.Dual, 0)).apply(allocator, star_b);
                    defer d_star_b.deinit(allocator);

                    const edge_volumes = state.mesh.simplices(1).items(.volume);
                    const dual_edge_volumes = state.mesh.dual_edge_volumes;
                    for (state.E.values, d_star_b.values, edge_volumes, dual_edge_volumes, state.J.values) |*electric, dsb, volume, dual_volume, current| {
                        std.debug.assert(dual_volume > 0.0);
                        electric.* += dt * ((volume / dual_volume) * dsb - current);
                    }
                },
                3 => {
                    var derivative = try (try state.operators.exteriorDerivative(cochain.Dual, 1)).apply(allocator, star_b);
                    defer derivative.deinit(allocator);

                    var curl_b = try (try state.operators.hodgeStarInverse(1)).apply(allocator, derivative);
                    defer curl_b.deinit(allocator);

                    for (state.E.values, curl_b.values, state.J.values) |*electric, curl_value, current| {
                        electric.* += dt * (curl_value - current);
                    }
                },
                else => unreachable,
            }
        }

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            return .{
                .state = try State.init(allocator, mesh),
            };
        }

        pub fn dipole(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: DipoleOptions(dim),
        ) !Self {
            comptime {
                if (dim != 2) @compileError("maxwell dipole is only defined for 2D");
            }

            var system = try Self.init(allocator, mesh);
            system.source = PointDipole(MeshType).init(mesh, options.frequency_hz, options.amplitude, options.center);
            system.boundary = options.boundary;
            return system;
        }

        pub fn cavity(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            options: CavityOptions(dim),
        ) !Self {
            var system = try Self.init(allocator, mesh);
            try initializeCavityState(dim, allocator, &system.state, options.time_step, options.extents);
            system.boundary = options.boundary;
            if (dim == 2) system.source = null;
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.state.deinit(allocator);
        }

        pub fn measurement(self: *const Self, allocator: std.mem.Allocator, which: Measurement) !f64 {
            return switch (which) {
                .energy => try electromagneticEnergy(allocator, &self.state),
            };
        }

        pub fn writeFields(self: *const Self, allocator: std.mem.Allocator, writer: anytype, fields: []const Field) !void {
            if (dim == 2) {
                var cell_data: [3]flux.io.DataArraySlice = undefined;
                var projected: [3]?[]f64 = .{ null, null, null };
                defer for (projected) |maybe_values| {
                    if (maybe_values) |values| allocator.free(values);
                };

                var count: usize = 0;
                for (fields) |field| {
                    cell_data[count] = try self.projectField2D(allocator, field, &projected[count]);
                    count += 1;
                }

                try flux.io.write(writer, self.state.mesh.*, &.{}, cell_data[0..count]);
                return;
            }

            var projected_fields: [3]common.viz.TetProjectionField = undefined;
            var count: usize = 0;
            for (fields) |field| {
                projected_fields[count] = self.projectField3D(field);
                count += 1;
            }

            try writeProjectedTetFields(count, allocator, writer, self.state.mesh, projected_fields);
        }

        fn projectField2D(
            self: *const Self,
            allocator: std.mem.Allocator,
            field: Field,
            projected: *?[]f64,
        ) !flux.io.DataArraySlice {
            return switch (field) {
                .electric => blk: {
                    projected.* = try flux.io.project_edges_to_faces(allocator, self.state.mesh.*, self.state.E.values);
                    break :blk .{ .name = "electric", .values = projected.*.? };
                },
                .magnetic => .{ .name = "magnetic", .values = self.state.B.values },
                .current => blk: {
                    projected.* = try flux.io.project_edges_to_faces(allocator, self.state.mesh.*, self.state.J.values);
                    break :blk .{ .name = "current", .values = projected.*.? };
                },
            };
        }

        fn projectField3D(self: *const Self, field: Field) common.viz.TetProjectionField {
            return switch (field) {
                .electric => .{ .name = "electric", .kind = .edge_abs_mean, .values = self.state.E.values },
                .magnetic => .{ .name = "magnetic", .kind = .face_abs_mean, .values = self.state.B.values },
                .current => .{ .name = "current", .kind = .edge_abs_mean, .values = self.state.J.values },
            };
        }
    };
}

fn primeOperators(
    comptime MeshType: type,
    operators: *operator_context_mod.OperatorContext(MeshType),
) !void {
    _ = try operators.exteriorDerivative(cochain.Primal, 1);
    _ = try operators.hodgeStar(2);

    switch (MeshType.topological_dimension) {
        2 => {
            _ = try operators.exteriorDerivative(cochain.Dual, 0);
            _ = try operators.hodgeStar(1);
        },
        3 => {
            _ = try operators.exteriorDerivative(cochain.Primal, 2);
            _ = try operators.exteriorDerivative(cochain.Dual, 1);
            _ = try operators.hodgeStarInverse(1);
        },
        else => unreachable,
    }
}

fn electromagneticEnergy(allocator: std.mem.Allocator, state: anytype) !f64 {
    const magnetic_degree = @TypeOf(state.B).degree;
    var star_e = try (try state.operators.hodgeStar(1)).apply(allocator, state.E);
    defer star_e.deinit(allocator);
    var star_b = try (try state.operators.hodgeStar(magnetic_degree)).apply(allocator, state.B);
    defer star_b.deinit(allocator);

    var electric_energy: f64 = 0.0;
    for (state.E.values, star_e.values) |lhs, rhs| {
        electric_energy += lhs * rhs;
    }

    var magnetic_energy: f64 = 0.0;
    for (state.B.values, star_b.values) |lhs, rhs| {
        magnetic_energy += lhs * rhs;
    }

    return 0.5 * (electric_energy + magnetic_energy);
}

fn PointDipole(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        edge_index: u32,
        edge_length: f64,
        frequency_hz: f64,
        amplitude: f64,

        pub fn init(mesh: *const MeshType, frequency_hz: f64, amplitude: f64, center: [MeshType.embedding_dimension]f64) Self {
            const simplex_1 = mesh.simplices(1);
            const edge_verts = simplex_1.items(.vertices);
            const lengths = simplex_1.items(.volume);
            const coords = mesh.vertices.slice().items(.coords);

            var best_edge: u32 = 0;
            var best_distance_squared = std.math.inf(f64);
            for (0..mesh.num_edges()) |edge_idx| {
                const v0 = coords[edge_verts[edge_idx][0]];
                const v1 = coords[edge_verts[edge_idx][1]];
                var midpoint: [MeshType.embedding_dimension]f64 = undefined;
                inline for (0..MeshType.embedding_dimension) |axis| {
                    midpoint[axis] = 0.5 * (v0[axis] + v1[axis]);
                }

                var distance_squared: f64 = 0.0;
                inline for (0..MeshType.embedding_dimension) |axis| {
                    const delta = midpoint[axis] - center[axis];
                    distance_squared += delta * delta;
                }
                if (distance_squared < best_distance_squared) {
                    best_distance_squared = distance_squared;
                    best_edge = @intCast(edge_idx);
                }
            }

            return .{
                .edge_index = best_edge,
                .edge_length = lengths[best_edge],
                .frequency_hz = frequency_hz,
                .amplitude = amplitude,
            };
        }

        pub fn apply(self: Self, current: anytype, time: f64) void {
            @memset(current.values, 0.0);
            current.values[self.edge_index] =
                (self.amplitude / self.edge_length) *
                @sin(2.0 * std.math.pi * self.frequency_hz * time);
        }
    };
}

pub fn DipoleOptions(comptime dim: u8) type {
    return struct {
        center: [dim]f64,
        frequency_hz: f64,
        amplitude: f64,
        boundary: BoundaryCondition = .pec,
    };
}

pub fn CavityOptions(comptime dim: u8) type {
    comptime {
        if (dim != 2 and dim != 3) @compileError("maxwell examples only support dimensions 2 and 3");
    }
    return struct {
        extents: [dim]f64,
        time_step: f64,
        boundary: BoundaryCondition = .pec,
    };
}

fn writeProjectedTetFields(
    field_count: usize,
    allocator: std.mem.Allocator,
    writer: anytype,
    mesh: anytype,
    projected_fields: [3]common.viz.TetProjectionField,
) !void {
    inline for (0..4) |count| {
        if (field_count == count) {
            if (count == 0) {
                return flux.io.write(writer, mesh.*, &.{}, &.{});
            }

            var fields: [count]common.viz.TetProjectionField = undefined;
            inline for (0..count) |field_idx| {
                fields[field_idx] = projected_fields[field_idx];
            }
            return common.viz.writeProjectedTetFields(count, allocator, writer, mesh, fields);
        }
    }
    unreachable;
}

pub fn CavityMeasurementProvider(comptime dim: u8) type {
    comptime {
        if (dim != 2 and dim != 3) @compileError("maxwell examples only support dimensions 2 and 3");
    }
    return struct {
        pub const Measurement = enum { energy, electric_l2, magnetic_l2 };

        extents: [dim]f64,
        time_step: f64,

        pub fn measurement(self: @This(), allocator: std.mem.Allocator, system: anytype, which: Measurement, time: f64) !f64 {
            return switch (which) {
                .energy => try system.measurement(allocator, .energy),
                .electric_l2 => try self.electricL2(allocator, system, time),
                .magnetic_l2 => try self.magneticL2(allocator, system, time),
            };
        }

        fn electricL2(self: @This(), allocator: std.mem.Allocator, system: anytype, time: f64) !f64 {
            var exact = try @TypeOf(system.state.E).init(allocator, system.state.mesh);
            defer exact.deinit(allocator);
            switch (dim) {
                2 => {
                    const Space = flux.forms.feec.WhitneySpace(@TypeOf(exact).MeshT, 1);
                    const project = bridges.DeRhamProjection(Space).init(Space.init(exact.mesh));
                    project.projectInto(&exact, TeMode(2, .electric, .{ 1, 0 }){
                        .time = time,
                        .extents = self.extents,
                    });
                },
                3 => {
                    const Space = flux.forms.feec.WhitneySpace(@TypeOf(exact).MeshT, 1);
                    const project = bridges.DeRhamProjection(Space).init(Space.init(exact.mesh));
                    project.projectInto(&exact, TmMode(3, .electric, .{ 1, 1, 0 }){
                        .time = time,
                        .extents = self.extents,
                    });
                },
                else => unreachable,
            }
            return relativeFormError(allocator, system.state.E, exact, system.state.operators, 1);
        }

        fn magneticL2(self: @This(), allocator: std.mem.Allocator, system: anytype, time: f64) !f64 {
            var exact = try @TypeOf(system.state.B).init(allocator, system.state.mesh);
            defer exact.deinit(allocator);
            switch (dim) {
                2 => {
                    const Space = flux.forms.feec.WhitneySpace(@TypeOf(exact).MeshT, 2);
                    const project = bridges.DeRhamProjection(Space).init(Space.init(exact.mesh));
                    project.projectInto(&exact, TeMode(2, .magnetic, .{ 1, 0 }){
                        .time = time - 0.5 * self.time_step,
                        .extents = self.extents,
                    });
                },
                3 => {
                    const Space = flux.forms.feec.WhitneySpace(@TypeOf(exact).MeshT, 2);
                    const project = bridges.DeRhamProjection(Space).init(Space.init(exact.mesh));
                    project.projectInto(&exact, TmMode(3, .magnetic, .{ 1, 1, 0 }){
                        .time = time - 0.5 * self.time_step,
                        .extents = self.extents,
                    });
                },
                else => unreachable,
            }
            return relativeFormError(allocator, system.state.B, exact, system.state.operators, 2);
        }
    };
}

fn relativeFormError(
    allocator: std.mem.Allocator,
    approx: anytype,
    exact: anytype,
    operator_context: anytype,
    comptime degree: u8,
) !f64 {
    var error_form = try @TypeOf(approx).init(allocator, approx.mesh);
    defer error_form.deinit(allocator);
    @memcpy(error_form.values, approx.values);
    error_form.sub(exact);

    var star_error = try (try operator_context.hodgeStar(degree)).apply(allocator, error_form);
    defer star_error.deinit(allocator);
    var star_exact = try (try operator_context.hodgeStar(degree)).apply(allocator, exact);
    defer star_exact.deinit(allocator);

    var error_energy: f64 = 0.0;
    for (error_form.values, star_error.values) |lhs, rhs| {
        error_energy += lhs * rhs;
    }

    var exact_energy: f64 = 0.0;
    for (exact.values, star_exact.values) |lhs, rhs| {
        exact_energy += lhs * rhs;
    }

    std.debug.assert(exact_energy > 0.0);
    return std.math.sqrt(error_energy / exact_energy);
}

const ModeQuantity = enum {
    electric,
    magnetic,
    potential,
};

fn modeWaveNumbers(comptime dim: u8, comptime indices: [dim]u32, extents: [dim]f64) [dim]f64 {
    var wave_numbers: [dim]f64 = undefined;
    inline for (0..dim) |axis| {
        wave_numbers[axis] = std.math.pi * @as(f64, @floatFromInt(indices[axis])) / extents[axis];
    }
    return wave_numbers;
}

fn modeAngularFrequency(comptime dim: u8, wave_numbers: [dim]f64) f64 {
    var sum_squares: f64 = 0.0;
    inline for (0..dim) |axis| {
        sum_squares += wave_numbers[axis] * wave_numbers[axis];
    }
    return std.math.sqrt(sum_squares);
}

fn TeMode(comptime dim: u8, comptime quantity: ModeQuantity, comptime indices: [dim]u32) type {
    comptime {
        if (dim != 2) @compileError("TE cavity modes are currently implemented only for 2D");
        if (indices[1] != 0) @compileError("2D TE cavity modes currently require indices = {m, 0}");
        if (quantity == .potential) @compileError("2D TE cavity mode does not expose a potential form in this example");
    }

    return struct {
        time: f64,
        extents: [dim]f64,

        pub fn evaluate(self: @This(), point: [dim]f64) if (quantity == .electric) [dim]f64 else f64 {
            const wave_numbers = modeWaveNumbers(dim, indices, self.extents);
            const omega = modeAngularFrequency(dim, wave_numbers);
            return switch (quantity) {
                .electric => .{
                    0.0,
                    @sin(wave_numbers[0] * point[0]) * @sin(omega * self.time),
                },
                .magnetic => @cos(wave_numbers[0] * point[0]) * @cos(omega * self.time),
                .potential => unreachable,
            };
        }
    };
}

fn TmMode(comptime dim: u8, comptime quantity: ModeQuantity, comptime indices: [dim]u32) type {
    comptime {
        if (dim != 3) @compileError("TM cavity modes are currently implemented only for 3D");
        if (indices[2] != 0) @compileError("3D TM cavity modes currently require indices = {m, n, 0}");
    }

    return struct {
        time: f64,
        extents: [dim]f64,

        pub fn evaluate(self: @This(), point: [dim]f64) if (quantity == .magnetic) [dim]f64 else [dim]f64 {
            const wave_numbers = modeWaveNumbers(dim, indices, self.extents);
            const omega = modeAngularFrequency(dim, wave_numbers);
            return switch (quantity) {
                .electric => .{
                    0.0,
                    0.0,
                    @sin(wave_numbers[0] * point[0]) * @sin(wave_numbers[1] * point[1]) * @sin(omega * self.time),
                },
                .potential => .{
                    0.0,
                    0.0,
                    (1.0 / omega) * @sin(wave_numbers[0] * point[0]) * @sin(wave_numbers[1] * point[1]) * @cos(omega * self.time),
                },
                .magnetic => .{
                    (wave_numbers[1] / omega) * @sin(wave_numbers[0] * point[0]) * @cos(wave_numbers[1] * point[1]) * @cos(omega * self.time),
                    -(wave_numbers[0] / omega) * @cos(wave_numbers[0] * point[0]) * @sin(wave_numbers[1] * point[1]) * @cos(omega * self.time),
                    0.0,
                },
            };
        }
    };
}

fn initializeCavityState(
    comptime dim: u8,
    allocator: std.mem.Allocator,
    state: anytype,
    time_step: f64,
    extents: [dim]f64,
) !void {
    if (dim == 2) {
        const Space = flux.forms.feec.WhitneySpace(@TypeOf(state.B).MeshT, 2);
        const project = bridges.DeRhamProjection(Space).init(Space.init(state.B.mesh));
        project.projectInto(&state.B, TeMode(2, .magnetic, .{ 1, 0 }){
            .time = -time_step / 2.0,
            .extents = extents,
        });
        return;
    }

    var potential = try @TypeOf(state.*).OneForm.init(allocator, state.mesh);
    defer potential.deinit(allocator);
    @memset(state.E.values, 0.0);
    {
        const Space = flux.forms.feec.WhitneySpace(@TypeOf(potential).MeshT, 1);
        const project = bridges.DeRhamProjection(Space).init(Space.init(potential.mesh));
        project.projectInto(&potential, TmMode(3, .potential, .{ 1, 1, 0 }){
            .time = -time_step / 2.0,
            .extents = extents,
        });
    }
    var exact_flux = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, potential);
    defer exact_flux.deinit(allocator);
    @memcpy(state.B.values, exact_flux.values);
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
