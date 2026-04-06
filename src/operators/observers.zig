//! Composable scalar diagnostics for simulation state.
//!
//! Observers are small values with a stable output name and an `evaluate`
//! method. Integrators or runners can evaluate any tuple of observers after a
//! step without knowing the physics-specific details of each diagnostic.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const exterior_derivative = @import("exterior_derivative.zig");
const wedge_product = @import("wedge_product.zig");

pub const Observation = struct {
    name: []const u8,
    value: f64,
};

pub const SignedCell = struct {
    cell_index: usize,
    orientation: i8,
};

fn validateStatePointer(comptime StatePointer: type) type {
    const info = @typeInfo(StatePointer);
    if (info != .pointer) {
        @compileError("observer evaluation requires a state pointer");
    }
    if (info.pointer.size != .one) {
        @compileError("observer evaluation requires a single-item state pointer");
    }
    return info.pointer.child;
}

fn validateObserverType(comptime ObserverType: type) void {
    if (!@hasDecl(ObserverType, "State")) {
        @compileError("observer type must declare `pub const State`");
    }
    if (!@hasDecl(ObserverType, "evaluate")) {
        @compileError("observer type must declare `pub fn evaluate(self, allocator, state, step)`");
    }
    if (!@hasField(ObserverType, "name")) {
        @compileError("observer type must have a `name` field");
    }
}

pub fn evaluateAll(
    allocator: std.mem.Allocator,
    state: anytype,
    step: usize,
    observers: anytype,
) ![observers.len]Observation {
    const StatePointer = @TypeOf(state);
    const State = comptime validateStatePointer(StatePointer);

    comptime {
        for (observers) |observer| {
            const ObserverType = @TypeOf(observer);
            validateObserverType(ObserverType);
            if (ObserverType.State != State) {
                @compileError("all observers in evaluateAll must share the same State type");
            }
        }
    }

    var results: [observers.len]Observation = undefined;
    inline for (observers, 0..) |observer, index| {
        results[index] = .{
            .name = observer.name,
            .value = try observer.evaluate(allocator, state, step),
        };
    }
    return results;
}

pub fn EnergyObserver(
    comptime StateType: type,
    comptime energyFn: fn (std.mem.Allocator, *const StateType, usize) anyerror!f64,
) type {
    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,

        pub fn evaluate(self: Self, allocator: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = self;
            return energyFn(allocator, state, step);
        }
    };
}

fn validateCochainType(comptime CochainType: type) void {
    if (!@hasDecl(CochainType, "MeshT") or !@hasDecl(CochainType, "degree") or !@hasDecl(CochainType, "duality")) {
        @compileError("observer requires a cochain type");
    }
}

pub fn L2NormObserver(
    comptime StateType: type,
    comptime CochainType: type,
    comptime fieldFn: fn (*const StateType) *const CochainType,
) type {
    comptime validateCochainType(CochainType);

    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,

        pub fn evaluate(self: Self, _: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = self;
            _ = step;
            const field = fieldFn(state);
            return std.math.sqrt(field.norm_squared());
        }
    };
}

pub fn MaxNormObserver(
    comptime StateType: type,
    comptime CochainType: type,
    comptime fieldFn: fn (*const StateType) *const CochainType,
) type {
    comptime validateCochainType(CochainType);

    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,

        pub fn evaluate(self: Self, _: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = self;
            _ = step;
            const field = fieldFn(state);
            var max_abs: f64 = 0.0;
            for (field.values) |value| {
                max_abs = @max(max_abs, @abs(value));
            }
            return max_abs;
        }
    };
}

pub fn CirculationObserver(
    comptime StateType: type,
    comptime CochainType: type,
    comptime fieldFn: fn (*const StateType) *const CochainType,
) type {
    comptime {
        validateCochainType(CochainType);
        if (CochainType.degree != 1) {
            @compileError("CirculationObserver requires a 1-cochain");
        }
    }

    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,
        loop: []const SignedCell,

        pub fn evaluate(self: Self, _: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = step;
            const field = fieldFn(state);
            var circulation: f64 = 0.0;
            for (self.loop) |term| {
                std.debug.assert(term.cell_index < field.values.len);
                std.debug.assert(term.orientation == -1 or term.orientation == 1);
                circulation += @as(f64, @floatFromInt(term.orientation)) * field.values[term.cell_index];
            }
            return circulation;
        }
    };
}

pub fn DivergenceNormObserver(
    comptime StateType: type,
    comptime CochainType: type,
    comptime fieldFn: fn (*const StateType) *const CochainType,
) type {
    comptime {
        validateCochainType(CochainType);
        if (CochainType.duality != cochain.Primal) {
            @compileError("DivergenceNormObserver currently requires a primal cochain");
        }
        if (CochainType.degree != CochainType.MeshT.topological_dimension - 1) {
            @compileError("DivergenceNormObserver requires an (n-1)-cochain so d(field) is an n-cochain");
        }
    }

    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,

        pub fn evaluate(self: Self, allocator: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = self;
            _ = step;
            const field = fieldFn(state);
            var divergence = try exterior_derivative.exterior_derivative(allocator, field.*);
            defer divergence.deinit(allocator);
            return std.math.sqrt(divergence.norm_squared());
        }
    };
}

pub fn HelicityObserver(
    comptime StateType: type,
    comptime CochainType: type,
    comptime fieldFn: fn (*const StateType) *const CochainType,
) type {
    comptime {
        validateCochainType(CochainType);
        if (CochainType.duality != cochain.Primal) {
            @compileError("HelicityObserver currently requires a primal cochain");
        }
        if (CochainType.MeshT.topological_dimension != 3) {
            @compileError("HelicityObserver is only defined on 3D meshes");
        }
        if (CochainType.degree != 1) {
            @compileError("HelicityObserver requires a primal 1-cochain");
        }
    }

    return struct {
        const Self = @This();
        pub const State = StateType;

        name: []const u8,

        pub fn evaluate(self: Self, allocator: std.mem.Allocator, state: *const StateType, step: usize) !f64 {
            _ = self;
            _ = step;
            const field = fieldFn(state);

            var derivative = try exterior_derivative.exterior_derivative(allocator, field.*);
            defer derivative.deinit(allocator);

            var helicity_density = try wedge_product.wedge(allocator, field.*, derivative);
            defer helicity_density.deinit(allocator);

            var helicity: f64 = 0.0;
            for (helicity_density.values) |value| {
                helicity += value;
            }
            return helicity;
        }
    };
}

const Mesh2D = topology.Mesh(2, 2);
const Mesh3D = topology.Mesh(3, 3);
const Primal0_2D = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const Primal1_2D = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const Primal2_2D = cochain.Cochain(Mesh2D, 2, cochain.Primal);
const Primal1_3D = cochain.Cochain(Mesh3D, 1, cochain.Primal);

const MockState2D = struct {
    scalar: Primal0_2D,
    flux: Primal1_2D,
};

const MockState3D = struct {
    velocity: Primal1_3D,
};

fn mockEnergy(_: std.mem.Allocator, state: *const MockState2D, step: usize) !f64 {
    return state.scalar.norm_squared() + state.flux.norm_squared() + @as(f64, @floatFromInt(step));
}

fn selectScalar(state: *const MockState2D) *const Primal0_2D {
    return &state.scalar;
}

fn selectFlux(state: *const MockState2D) *const Primal1_2D {
    return &state.flux;
}

fn selectVelocity(state: *const MockState3D) *const Primal1_3D {
    return &state.velocity;
}

test "evaluateAll composes energy, L2, and max observers on a mock state" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var scalar = try Primal0_2D.init(allocator, &mesh);
    defer scalar.deinit(allocator);
    var flux = try Primal1_2D.init(allocator, &mesh);
    defer flux.deinit(allocator);

    for (scalar.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index + 1));
    }
    for (flux.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index)) - 2.5;
    }

    const state = MockState2D{
        .scalar = scalar,
        .flux = flux,
    };

    const Energy = EnergyObserver(MockState2D, mockEnergy);
    const L2 = L2NormObserver(MockState2D, Primal1_2D, selectFlux);
    const Max = MaxNormObserver(MockState2D, Primal1_2D, selectFlux);

    const results = try evaluateAll(allocator, &state, 7, .{
        Energy{ .name = "energy" },
        L2{ .name = "flux_l2" },
        Max{ .name = "flux_max" },
    });

    try testing.expectEqualStrings("energy", results[0].name);
    try testing.expectEqualStrings("flux_l2", results[1].name);
    try testing.expectEqualStrings("flux_max", results[2].name);
    try testing.expectApproxEqAbs(try mockEnergy(allocator, &state, 7), results[0].value, 1e-12);
    try testing.expectApproxEqAbs(std.math.sqrt(state.flux.norm_squared()), results[1].value, 1e-12);

    var expected_max: f64 = 0.0;
    for (state.flux.values) |value| {
        expected_max = @max(expected_max, @abs(value));
    }
    try testing.expectApproxEqAbs(expected_max, results[2].value, 1e-12);
}

test "circulation observer sums an oriented loop through a 1-cochain" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var scalar = try Primal0_2D.init(allocator, &mesh);
    defer scalar.deinit(allocator);
    var flux = try Primal1_2D.init(allocator, &mesh);
    defer flux.deinit(allocator);
    for (flux.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index + 1));
    }

    const state = MockState2D{
        .scalar = scalar,
        .flux = flux,
    };
    const loop = [_]SignedCell{
        .{ .cell_index = 0, .orientation = 1 },
        .{ .cell_index = 1, .orientation = -1 },
        .{ .cell_index = 4, .orientation = 1 },
    };

    const Circulation = CirculationObserver(MockState2D, Primal1_2D, selectFlux);
    const observer = Circulation{
        .name = "circulation",
        .loop = &loop,
    };

    const value = try observer.evaluate(allocator, &state, 0);
    try testing.expectApproxEqAbs(flux.values[0] - flux.values[1] + flux.values[4], value, 1e-12);
}

test "divergence norm observer vanishes on exact top-degree derivative by dd = 0" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var potential = try Primal0_2D.init(allocator, &mesh);
    defer potential.deinit(allocator);
    for (potential.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index * 3 + 1));
    }

    var exact_flux = try exterior_derivative.exterior_derivative(allocator, potential);
    defer exact_flux.deinit(allocator);

    var scalar = try Primal0_2D.init(allocator, &mesh);
    defer scalar.deinit(allocator);
    const state = MockState2D{
        .scalar = scalar,
        .flux = exact_flux,
    };

    const Divergence = DivergenceNormObserver(MockState2D, Primal1_2D, selectFlux);
    const observer = Divergence{ .name = "divergence_l2" };

    const value = try observer.evaluate(allocator, &state, 0);
    try testing.expectApproxEqAbs(@as(f64, 0.0), value, 1e-12);
}

test "helicity observer matches manual sum of u wedge du on tetrahedral mesh" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var velocity = try Primal1_3D.init(allocator, &mesh);
    defer velocity.deinit(allocator);
    for (velocity.values, 0..) |*value, index| {
        value.* = 0.25 * @as(f64, @floatFromInt(index + 1));
    }

    const state = MockState3D{ .velocity = velocity };

    const Helicity = HelicityObserver(MockState3D, Primal1_3D, selectVelocity);
    const observer = Helicity{ .name = "helicity" };

    const observed = try observer.evaluate(allocator, &state, 0);

    var derivative = try exterior_derivative.exterior_derivative(allocator, state.velocity);
    defer derivative.deinit(allocator);
    var density = try wedge_product.wedge(allocator, state.velocity, derivative);
    defer density.deinit(allocator);

    var expected: f64 = 0.0;
    for (density.values) |value| {
        expected += value;
    }

    try testing.expectApproxEqAbs(expected, observed, 1e-12);
}

test "divergence observer matches manual exterior-derivative norm" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var scalar = try Primal0_2D.init(allocator, &mesh);
    defer scalar.deinit(allocator);
    var flux = try Primal1_2D.init(allocator, &mesh);
    defer flux.deinit(allocator);
    for (flux.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index + 2));
    }

    const state = MockState2D{
        .scalar = scalar,
        .flux = flux,
    };

    const Divergence = DivergenceNormObserver(MockState2D, Primal1_2D, selectFlux);
    const observer = Divergence{ .name = "divergence_l2" };
    const observed = try observer.evaluate(allocator, &state, 0);

    var divergence = try exterior_derivative.exterior_derivative(allocator, state.flux);
    defer divergence.deinit(allocator);
    const expected = std.math.sqrt(divergence.norm_squared());
    try testing.expectApproxEqAbs(expected, observed, 1e-12);
}

test "manual wedge for helicity produces a top form" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var velocity = try Primal1_3D.init(allocator, &mesh);
    defer velocity.deinit(allocator);
    for (velocity.values, 0..) |*value, index| {
        value.* = @as(f64, @floatFromInt(index + 1));
    }

    var derivative = try exterior_derivative.exterior_derivative(allocator, velocity);
    defer derivative.deinit(allocator);
    var density = try wedge_product.wedge(allocator, velocity, derivative);
    defer density.deinit(allocator);

    try testing.expectEqual(@as(usize, mesh.num_tets()), density.values.len);
    _ = Primal2_2D;
}
