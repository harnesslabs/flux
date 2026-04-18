const std = @import("std");
const flux = @import("flux");

const bridges = flux.operators.bridges;
const cochain = flux.forms;
const feec_context_mod = flux.operators.feec.context;
const linear_system = flux.math.linear_system;
const sparse = flux.math.sparse;

pub const SurfaceKind = enum {
    plane,
    sphere,
};

const InitialCondition = enum {
    sine_mode,
    z_mode,
};

pub fn Diffusion(comptime surface_kind: SurfaceKind, comptime MeshType: type) type {
    return struct {
        const Self = @This();

        pub const Temperature = cochain.Cochain(MeshType, 0, cochain.Primal);
        pub const Field = enum { temperature };

        mesh: *const MeshType,
        runtime: if (surface_kind == .plane) PlaneRuntime(MeshType) else SphereRuntime(MeshType),
        temperature: Temperature,

        pub const BackwardEuler = struct {
            pub fn initialize(_: std.mem.Allocator, system: *Self, _: f64) void {
                system.runtime.linear_system.seedSolution(system.temperature.values);
            }

            pub fn advance(allocator: std.mem.Allocator, system: *Self, _: f64) !void {
                _ = allocator;
                try system.stepBackwardEuler();
            }
        };

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            time_step: f64,
        ) !Self {
            var runtime = if (surface_kind == .plane)
                try PlaneRuntime(MeshType).init(allocator, mesh, time_step)
            else
                try SphereRuntime(MeshType).init(allocator, mesh, time_step);
            errdefer runtime.deinit(allocator);

            var temperature = try Temperature.init(allocator, mesh);
            errdefer temperature.deinit(allocator);

            return .{
                .mesh = mesh,
                .runtime = runtime,
                .temperature = temperature,
            };
        }

        pub fn plane(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            time_step: f64,
        ) !Self {
            comptime {
                if (surface_kind != .plane) @compileError("plane constructor only applies to plane diffusion");
            }
            var system = try Self.init(allocator, mesh, time_step);
            projectScalarInto(&system.temperature, PlaneExactMode{ .time = 0.0 });
            return system;
        }

        pub fn sphere(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            time_step: f64,
        ) !Self {
            comptime {
                if (surface_kind != .sphere) @compileError("sphere constructor only applies to surface diffusion");
            }
            var system = try Self.init(allocator, mesh, time_step);
            projectScalarInto(&system.temperature, SphereExactMode{ .time = 0.0 });
            return system;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.temperature.deinit(allocator);
            self.runtime.deinit(allocator);
        }

        pub fn writeFields(self: *const Self, allocator: std.mem.Allocator, writer: anytype, fields: []const Field) !void {
            _ = allocator;
            std.debug.assert(fields.len <= 1);
            if (fields.len == 0) {
                return flux.io.write(writer, self.mesh.*, &.{}, &.{});
            }

            try flux.io.write(writer, self.mesh.*, &.{
                .{ .name = "temperature", .values = self.temperature.values },
            }, &.{});
        }

        fn stepBackwardEuler(self: *Self) !void {
            if (surface_kind == .plane) {
                const masses = self.mesh.vertices.slice().items(.dual_volume);
                const full_rhs = self.runtime.linear_system.fullRhsValues();
                std.debug.assert(full_rhs.len == self.temperature.values.len);

                for (full_rhs, masses, self.temperature.values) |*rhs_value, mass, state_value| {
                    rhs_value.* = mass * state_value;
                }

                _ = try self.runtime.linear_system.solveHomogeneous(self.temperature.values);
                return;
            }

            const rhs = self.runtime.linear_system.rhsValues();
            const solution = self.runtime.linear_system.solutionValues();
            std.debug.assert(rhs.len == self.temperature.values.len);
            std.debug.assert(solution.len == self.temperature.values.len);

            for (rhs, self.temperature.values, self.runtime.masses) |*rhs_value, state_value, mass| {
                rhs_value.* = mass * state_value;
            }

            _ = try self.runtime.linear_system.solve();
            @memcpy(self.temperature.values, solution);
        }
    };
}

fn PlaneRuntime(comptime MeshType: type) type {
    return struct {
        operator_context: *feec_context_mod.OperatorContext(MeshType),
        linear_system: linear_system.LinearSystem,

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            time_step: f64,
        ) !@This() {
            const operator_context = try feec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
            errdefer operator_context.deinit();

            const laplacian = try operator_context.laplacian(0);
            const stiffness = laplacian.stiffness;
            const masses = mesh.vertices.slice().items(.dual_volume);
            const elimination_map = try linear_system.EliminationMap.initBoundary(MeshType, allocator, mesh, 0);

            var triplets = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
            defer triplets.deinit(allocator);
            for (0..mesh.num_vertices()) |row_idx_usize| {
                const row_idx: u32 = @intCast(row_idx_usize);
                try triplets.addEntry(allocator, row_idx, row_idx, masses[row_idx]);

                const row = stiffness.row(row_idx);
                for (row.cols, row.vals) |col_idx, value| {
                    try triplets.addEntry(allocator, row_idx, col_idx, time_step * value);
                }
            }

            var full_matrix = try triplets.build(allocator);
            errdefer full_matrix.deinit(allocator);

            var system_runtime = try linear_system.LinearSystem.eliminate(
                allocator,
                full_matrix,
                elimination_map,
                .{},
            );
            errdefer system_runtime.deinit(allocator);
            full_matrix.deinit(allocator);

            return .{
                .operator_context = operator_context,
                .linear_system = system_runtime,
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.linear_system.deinit(allocator);
            self.operator_context.deinit();
        }
    };
}

fn SphereRuntime(comptime MeshType: type) type {
    return struct {
        operator_context: *feec_context_mod.OperatorContext(MeshType),
        linear_system: linear_system.LinearSystem,
        masses: []f64,

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            time_step: f64,
        ) !@This() {
            const operator_context = try feec_context_mod.OperatorContext(MeshType).init(allocator, mesh);
            errdefer operator_context.deinit();
            const stiffness = (try operator_context.laplacian(0)).stiffness;

            const masses = try assembleLumpedSurfaceMasses(allocator, mesh);
            errdefer allocator.free(masses);

            var assembler = sparse.TripletAssembler(f64).init(mesh.num_vertices(), mesh.num_vertices());
            defer assembler.deinit(allocator);

            for (masses, 0..) |mass, row_idx| {
                try assembler.addEntry(allocator, @intCast(row_idx), @intCast(row_idx), mass);
                const row = stiffness.row(@intCast(row_idx));
                for (row.cols, row.vals) |col_idx, value| {
                    try assembler.addEntry(allocator, @intCast(row_idx), col_idx, time_step * value);
                }
            }

            var system_matrix = try assembler.build(allocator);
            errdefer system_matrix.deinit(allocator);

            var system_runtime = try linear_system.LinearSystem.init(
                allocator,
                system_matrix,
                .{},
            );
            errdefer system_runtime.deinit(allocator);

            return .{
                .operator_context = operator_context,
                .linear_system = system_runtime,
                .masses = masses,
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.linear_system.deinit(allocator);
            allocator.free(self.masses);
            self.operator_context.deinit();
        }
    };
}

pub fn ReferenceMeasurementProvider(comptime surface_kind: SurfaceKind, comptime MeshType: type) type {
    _ = MeshType;
    return struct {
        pub const Measurement = enum { l2_error };

        pub fn measurement(_: @This(), allocator: std.mem.Allocator, system: anytype, which: Measurement, time: f64) !f64 {
            _ = which;
            var exact = try @TypeOf(system.temperature).init(allocator, system.mesh);
            defer exact.deinit(allocator);

            if (surface_kind == .plane) {
                projectScalarInto(&exact, PlaneExactMode{ .time = time });
            } else {
                projectScalarInto(&exact, SphereExactMode{ .time = time });
            }
            return weightedL2Error(surface_kind, system.mesh, system.temperature.values, exact.values);
        }
    };
}

const PlaneExactMode = struct {
    time: f64,

    pub fn evaluate(self: @This(), point: [2]f64) f64 {
        const eigenvalue = 2.0 * std.math.pi * std.math.pi;
        return std.math.exp(-eigenvalue * self.time) *
            std.math.sin(std.math.pi * point[0]) *
            std.math.sin(std.math.pi * point[1]);
    }
};

const SphereExactMode = struct {
    time: f64,

    pub fn evaluate(self: @This(), point: [3]f64) f64 {
        return std.math.exp(-2.0 * self.time) * point[2];
    }
};

fn projectScalarInto(target: anytype, continuous_form: anytype) void {
    const CochainType = @TypeOf(target.*);
    const Space = flux.forms.feec.WhitneySpace(CochainType.MeshT, 0);
    const projector = bridges.DeRhamProjection(Space).init(Space.init(target.mesh));
    projector.projectInto(target, continuous_form);
}

fn weightedL2Error(
    comptime surface_kind: SurfaceKind,
    mesh: anytype,
    approx: []const f64,
    exact: []const f64,
) f64 {
    std.debug.assert(approx.len == exact.len);
    if (surface_kind == .plane) {
        const dual_volumes = mesh.vertices.slice().items(.dual_volume);
        var error_sq: f64 = 0.0;
        var measure: f64 = 0.0;
        for (approx, exact, dual_volumes) |approx_value, exact_value, dual_volume| {
            const diff = approx_value - exact_value;
            error_sq += diff * diff * dual_volume;
            measure += dual_volume;
        }
        return @sqrt(error_sq / measure);
    }

    const face_vertices = mesh.simplices(2).items(.vertices);
    const face_areas = mesh.simplices(2).items(.volume);
    var error_sq: f64 = 0.0;
    var measure: f64 = 0.0;
    for (face_vertices, face_areas) |face, area| {
        const lumped = area / 3.0;
        const diff0 = approx[face[0]] - exact[face[0]];
        const diff1 = approx[face[1]] - exact[face[1]];
        const diff2 = approx[face[2]] - exact[face[2]];
        error_sq += lumped * (diff0 * diff0 + diff1 * diff1 + diff2 * diff2);
        measure += 3.0 * lumped;
    }
    return @sqrt(error_sq / measure);
}

fn assembleLumpedSurfaceMasses(
    allocator: std.mem.Allocator,
    mesh: anytype,
) ![]f64 {
    const masses = try allocator.alloc(f64, mesh.num_vertices());
    @memset(masses, 0.0);

    const face_vertices = mesh.simplices(2).items(.vertices);
    const face_areas = mesh.simplices(2).items(.volume);
    for (face_vertices, face_areas) |face, area| {
        const lumped = area / 3.0;
        masses[face[0]] += lumped;
        masses[face[1]] += lumped;
        masses[face[2]] += lumped;
    }
    return masses;
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
