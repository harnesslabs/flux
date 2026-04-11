//! Library-side evolution mechanics for reusable state-stepping systems.
//!
//! This layer owns stepper lifetime and repeated timestep orchestration state
//! that sit between assembled operators and higher-level example run loops.
//! Exact/reference evolution is one consumer, not the only supported shape.

const std = @import("std");
const testing = std.testing;

fn validateStepperType(comptime StepperType: type) void {
    if (!@hasDecl(StepperType, "step")) {
        @compileError("evolution builder stepper must declare `pub fn step(self: *Stepper, allocator) !void`");
    }
    if (!@hasDecl(StepperType, "deinit")) {
        @compileError("evolution stepper must declare `pub fn deinit(self: *Stepper, allocator) void`");
    }

    const step_info = @typeInfo(@TypeOf(StepperType.step));
    if (step_info != .@"fn") {
        @compileError("evolution stepper `step` must be a function");
    }
    const step_fn_info = step_info.@"fn";
    if (step_fn_info.params.len != 2) {
        @compileError("evolution stepper `step` must take (self: *Stepper, allocator)");
    }
    if (step_fn_info.params[0].type != *StepperType) {
        @compileError("evolution stepper `step` receiver must be *Stepper");
    }
    if (step_fn_info.params[1].type != std.mem.Allocator) {
        @compileError("evolution stepper `step` allocator parameter must be std.mem.Allocator");
    }
}

fn validateBuilderType(comptime StateType: type, comptime BuilderType: type) void {
    if (!@hasDecl(BuilderType, "Stepper")) {
        @compileError("evolution builder must declare `pub const Stepper`");
    }
    if (!@hasDecl(BuilderType, "initStepper")) {
        @compileError("evolution builder must declare `pub fn initStepper(self, allocator, state) !Stepper`");
    }

    const StepperType = BuilderType.Stepper;
    validateStepperType(StepperType);

    const init_info = @typeInfo(@TypeOf(BuilderType.initStepper));
    if (init_info != .@"fn") {
        @compileError("evolution builder `initStepper` must be a function");
    }
    const init_fn_info = init_info.@"fn";
    if (init_fn_info.params.len != 3) {
        @compileError("evolution builder `initStepper` must take (self, allocator, state)");
    }
    if (init_fn_info.params[1].type != std.mem.Allocator) {
        @compileError("evolution builder `initStepper` allocator parameter must be std.mem.Allocator");
    }
    if (init_fn_info.params[2].type != StateType) {
        @compileError("evolution builder `initStepper` state parameter must match the evolution state type");
    }
}

fn validateExactInitializerType(comptime MeshType: type, comptime Value: type, comptime ExactInitializerType: type) void {
    if (!@hasDecl(ExactInitializerType, "fill")) {
        @compileError("exact initializer must declare `pub fn fill(self, mesh, values, time) void`");
    }

    const info = @typeInfo(@TypeOf(ExactInitializerType.fill));
    if (info != .@"fn") {
        @compileError("exact initializer `fill` must be a function");
    }
    const fn_info = info.@"fn";
    if (fn_info.params.len != 4) {
        @compileError("exact initializer `fill` must take (self, *const MeshType, []Value, f64)");
    }
    if (fn_info.params[1].type != *const MeshType) {
        @compileError("exact initializer `fill` parameter 1 must be *const MeshType");
    }
    if (fn_info.params[2].type != []Value) {
        @compileError("exact initializer `fill` parameter 2 must be []Value");
    }
    if (fn_info.params[3].type != f64) {
        @compileError("exact initializer `fill` parameter 3 must be f64");
    }
}

fn validateStepFunction(comptime StateType: type, comptime StepFunctionType: type) void {
    const info = @typeInfo(StepFunctionType);
    if (info != .@"fn") {
        @compileError("fixed-time evolution step function must be a function");
    }
    const fn_info = info.@"fn";
    if (fn_info.params.len != 3) {
        @compileError("fixed-time evolution step function must take (allocator, state, dt)");
    }
    if (fn_info.params[0].type != std.mem.Allocator) {
        @compileError("fixed-time evolution step function allocator parameter must be std.mem.Allocator");
    }
    if (fn_info.params[1].type != StateType) {
        @compileError("fixed-time evolution step function state parameter must match the evolution state type");
    }
    if (fn_info.params[2].type != f64) {
        @compileError("fixed-time evolution step function dt parameter must be f64");
    }
}

fn validateErrorMeasureType(comptime MeshType: type, comptime Value: type, comptime ErrorMeasureType: type) void {
    if (!@hasDecl(ErrorMeasureType, "compute")) {
        @compileError("error measure must declare `pub fn compute(self, mesh, approx, exact) f64`");
    }

    const info = @typeInfo(@TypeOf(ErrorMeasureType.compute));
    if (info != .@"fn") {
        @compileError("error measure `compute` must be a function");
    }
    const fn_info = info.@"fn";
    if (fn_info.params.len != 4) {
        @compileError("error measure `compute` must take (self, *const MeshType, []const Value, []const Value)");
    }
    if (fn_info.params[1].type != *const MeshType) {
        @compileError("error measure `compute` parameter 1 must be *const MeshType");
    }
    if (fn_info.params[2].type != []const Value) {
        @compileError("error measure `compute` parameter 2 must be []const Value");
    }
    if (fn_info.params[3].type != []const Value) {
        @compileError("error measure `compute` parameter 3 must be []const Value");
    }
}

/// Evolution object for reusable state stepping.
///
/// The caller owns the primary state storage. The evolution object owns the
/// stepper instance, including any scratch that stepper embeds.
pub fn Evolution(
    comptime StateType: type,
    comptime StepperType: type,
) type {
    comptime {
        validateStepperType(StepperType);
    }

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        state: StateType,
        stepper: StepperType,
        step_count: u32,

        pub fn init(
            allocator: std.mem.Allocator,
            state: StateType,
            stepper: StepperType,
        ) Self {
            return .{
                .allocator = allocator,
                .state = state,
                .stepper = stepper,
                .step_count = 0,
            };
        }

        pub fn deinit(self: *Self) void {
            self.stepper.deinit(self.allocator);
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try self.stepper.step(allocator);
            self.step_count += 1;
        }

        pub fn currentState(self: *const Self) StateType {
            return self.state;
        }

        pub fn stepCount(self: *const Self) u32 {
            return self.step_count;
        }
    };
}

/// Evolution object for strategies that only need `(state, dt)` and own no scratch.
pub fn FixedTimeEvolution(
    comptime StateType: type,
    comptime step_function: anytype,
) type {
    comptime {
        validateStepFunction(StateType, @TypeOf(step_function));
    }

    return struct {
        const Self = @This();

        state: StateType,
        dt: f64,
        step_count: u32,

        pub fn init(state: StateType, dt: f64) Self {
            return .{
                .state = state,
                .dt = dt,
                .step_count = 0,
            };
        }

        pub fn deinit(self: *Self) void {
            _ = self;
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try step_function(allocator, self.state, self.dt);
            self.step_count += 1;
        }

        pub fn currentState(self: *const Self) StateType {
            return self.state;
        }

        pub fn stepCount(self: *const Self) u32 {
            return self.step_count;
        }
    };
}

/// Evolution object for a single time-stepped scalar field with an exact/reference solution.
///
/// The caller owns the mesh and the primary state slice. The evolution object
/// owns the reusable exact/reference buffer and the builder-produced stepper.
pub fn ExactEvolution(
    comptime MeshType: type,
    comptime Value: type,
    comptime BuilderType: type,
    comptime ExactInitializerType: type,
    comptime ErrorMeasureType: type,
) type {
    comptime {
        validateBuilderType([]Value, BuilderType);
        validateExactInitializerType(MeshType, Value, ExactInitializerType);
        validateErrorMeasureType(MeshType, Value, ErrorMeasureType);
    }

    return struct {
        const Self = @This();
        const BaseEvolution = Evolution([]Value, BuilderType.Stepper);

        mesh: *const MeshType,
        base: BaseEvolution,
        exact_values: []Value,
        exact_initializer: ExactInitializerType,
        error_measure: ErrorMeasureType,

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            state_values: []Value,
            builder: BuilderType,
            exact_initializer: ExactInitializerType,
            error_measure: ErrorMeasureType,
        ) !Self {
            const exact_values = try allocator.alloc(Value, state_values.len);
            errdefer allocator.free(exact_values);

            const stepper = try builder.initStepper(allocator, state_values);
            errdefer {
                var stepper_local = stepper;
                stepper_local.deinit(allocator);
            }
            const base = BaseEvolution.init(allocator, state_values, stepper);

            return .{
                .mesh = mesh,
                .base = base,
                .exact_values = exact_values,
                .exact_initializer = exact_initializer,
                .error_measure = error_measure,
            };
        }

        pub fn deinit(self: *Self) void {
            self.base.deinit();
            self.base.allocator.free(self.exact_values);
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try self.base.step(allocator);
        }

        pub fn fillExact(self: *Self, time: f64) void {
            self.exact_initializer.fill(self.mesh, self.exact_values, time);
        }

        pub fn l2Error(self: *const Self) f64 {
            return self.error_measure.compute(self.mesh, self.base.currentState(), self.exact_values);
        }

        pub fn stateValues(self: *const Self) []const Value {
            return self.base.currentState();
        }

        pub fn exactValues(self: *Self) []Value {
            return self.exact_values;
        }

        pub fn stepCount(self: *const Self) u32 {
            return self.base.stepCount();
        }
    };
}

const MockMesh = struct {};

const MockStepper = struct {
    state_values: []f64,
    rhs: []f64,
    solution: []f64,

    pub fn step(self: *@This(), _: std.mem.Allocator) !void {
        for (self.state_values, self.rhs, self.solution) |*state_value, *rhs_value, *solution_value| {
            rhs_value.* = state_value.*;
            solution_value.* = rhs_value.* + 1.0;
            state_value.* = solution_value.*;
        }
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.solution);
        allocator.free(self.rhs);
    }
};

const MockBuilder = struct {
    pub const Stepper = MockStepper;

    scratch_len: usize,

    pub fn initStepper(self: @This(), allocator: std.mem.Allocator, state_values: []f64) !Stepper {
        const stepper = Stepper{
            .state_values = state_values,
            .rhs = try allocator.alloc(f64, self.scratch_len),
            .solution = try allocator.alloc(f64, self.scratch_len),
        };
        @memcpy(stepper.solution, state_values);
        return stepper;
    }
};

const MockExactInitializer = struct {
    pub fn fill(_: @This(), _: *const MockMesh, values: []f64, time: f64) void {
        for (values, 0..) |*value, idx| {
            value.* = time + @as(f64, @floatFromInt(idx));
        }
    }
};

const MockErrorMeasure = struct {
    pub fn compute(_: @This(), _: *const MockMesh, approx: []const f64, exact: []const f64) f64 {
        var sum: f64 = 0.0;
        for (approx, exact) |approx_value, exact_value| {
            sum += @abs(approx_value - exact_value);
        }
        return sum;
    }
};

test "Evolution owns direct stepper state" {
    const TestEvolution = Evolution([]f64, MockStepper);
    const allocator = testing.allocator;
    var state = [_]f64{ 1.0, 2.0, 3.0 };
    const builder = MockBuilder{ .scratch_len = state.len };
    const stepper = try builder.initStepper(allocator, state[0..]);
    var evolution = TestEvolution.init(allocator, state[0..], stepper);
    defer evolution.deinit();

    try testing.expectEqual(state.len, evolution.stepper.rhs.len);
    try testing.expectEqual(state.len, evolution.stepper.solution.len);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());

    try evolution.step(allocator);
    try testing.expectEqual(@as(u32, 1), evolution.stepCount());
    try testing.expectApproxEqAbs(2.0, state[0], 1e-15);
    try testing.expectApproxEqAbs(3.0, state[1], 1e-15);
    try testing.expectApproxEqAbs(4.0, state[2], 1e-15);
}

test "FixedTimeEvolution wraps a plain state-step function" {
    const StepFunction = struct {
        pub fn step(_: std.mem.Allocator, state: []f64, dt: f64) !void {
            for (state) |*value| value.* += dt;
        }
    }.step;

    const TestEvolution = FixedTimeEvolution([]f64, StepFunction);
    var state = [_]f64{ 1.0, 2.0 };
    var evolution = TestEvolution.init(state[0..], 0.5);
    defer evolution.deinit();

    try evolution.step(testing.allocator);
    try testing.expectEqual(@as(u32, 1), evolution.stepCount());
    try testing.expectApproxEqAbs(1.5, state[0], 1e-15);
    try testing.expectApproxEqAbs(2.5, state[1], 1e-15);
}

test "ExactEvolution owns exact buffer and builder-produced stepper state" {
    const TestEvolution = ExactEvolution(MockMesh, f64, MockBuilder, MockExactInitializer, MockErrorMeasure);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0, 3.0 };

    var evolution = try TestEvolution.init(
        allocator,
        &mesh,
        state[0..],
        .{ .scratch_len = state.len },
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer evolution.deinit();

    try testing.expectEqual(state.len, evolution.exactValues().len);
    try testing.expectEqual(state.len, evolution.base.stepper.rhs.len);
    try testing.expectEqual(state.len, evolution.base.stepper.solution.len);
}

test "ExactEvolution computes error against the current exact field" {
    const TestEvolution = ExactEvolution(MockMesh, f64, MockBuilder, MockExactInitializer, MockErrorMeasure);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };

    var evolution = try TestEvolution.init(
        allocator,
        &mesh,
        state[0..],
        .{ .scratch_len = state.len },
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer evolution.deinit();

    evolution.fillExact(0.0);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(2.0, evolution.l2Error(), 1e-15);
}
