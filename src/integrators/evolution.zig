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

/// Evolution object for reusable state stepping.
///
/// The caller owns the primary state storage. The evolution object owns the
/// stepper instance, including any scratch that stepper embeds.
pub fn Evolution(
    comptime StateType: type,
    comptime StepperType: type,
    comptime AuxType: type,
) type {
    comptime {
        validateStepperType(StepperType);
    }

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        state: StateType,
        stepper: StepperType,
        aux: AuxType,
        step_count: u32,

        pub fn init(
            allocator: std.mem.Allocator,
            state: StateType,
            stepper: StepperType,
            aux: AuxType,
        ) Self {
            return .{
                .allocator = allocator,
                .state = state,
                .stepper = stepper,
                .aux = aux,
                .step_count = 0,
            };
        }

        pub fn deinit(self: *Self) void {
            self.stepper.deinit(self.allocator);
            if (AuxType != void) {
                self.aux.deinit(self.allocator);
            }
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try self.stepper.step(allocator);
            self.step_count += 1;
        }

        pub fn currentState(self: *const Self) StateType {
            return self.state;
        }

        pub fn stateValues(self: *const Self) StateType {
            return self.state;
        }

        pub fn fillExact(self: *Self, time: f64) void {
            if (!@hasDecl(AuxType, "fillExact")) {
                @compileError("Evolution auxiliary state must declare `pub fn fillExact(self, time) void` to support exact/reference runs");
            }
            self.aux.fillExact(time);
        }

        pub fn exactValues(self: *Self) []f64 {
            if (!@hasDecl(AuxType, "exactValues")) {
                @compileError("Evolution auxiliary state must declare `pub fn exactValues(self)` to expose reference values");
            }
            return self.aux.exactValues();
        }

        pub fn l2Error(self: *const Self) f64 {
            if (!@hasDecl(AuxType, "l2Error")) {
                @compileError("Evolution auxiliary state must declare `pub fn l2Error(self, state)` to evaluate errors");
            }
            return self.aux.l2Error(self.state);
        }

        pub fn stepCount(self: *const Self) u32 {
            return self.step_count;
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

const MockExactAux = struct {
    mesh: *const MockMesh,
    exact_values: []f64,
    initializer: MockExactInitializer,
    error_measure: MockErrorMeasure,

    pub fn init(
        allocator: std.mem.Allocator,
        mesh: *const MockMesh,
        len: usize,
        initializer: MockExactInitializer,
        error_measure: MockErrorMeasure,
    ) !@This() {
        return .{
            .mesh = mesh,
            .exact_values = try allocator.alloc(f64, len),
            .initializer = initializer,
            .error_measure = error_measure,
        };
    }

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.exact_values);
    }

    pub fn fillExact(self: *@This(), time: f64) void {
        self.initializer.fill(self.mesh, self.exact_values, time);
    }

    pub fn exactValues(self: *@This()) []f64 {
        return self.exact_values;
    }

    pub fn l2Error(self: *const @This(), state_values: []const f64) f64 {
        return self.error_measure.compute(self.mesh, state_values, self.exact_values);
    }
};

test "Evolution owns direct stepper state" {
    const TestEvolution = Evolution([]f64, MockStepper, void);
    const allocator = testing.allocator;
    var state = [_]f64{ 1.0, 2.0, 3.0 };
    const builder = MockBuilder{ .scratch_len = state.len };
    const stepper = try builder.initStepper(allocator, state[0..]);
    var evolution = TestEvolution.init(allocator, state[0..], stepper, {});
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

test "Evolution owns exact buffer through auxiliary state" {
    const TestEvolution = Evolution([]f64, MockStepper, MockExactAux);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0, 3.0 };
    const builder = MockBuilder{ .scratch_len = state.len };
    const stepper = try builder.initStepper(allocator, state[0..]);
    const aux = try MockExactAux.init(
        allocator,
        &mesh,
        state.len,
        MockExactInitializer{},
        MockErrorMeasure{},
    );

    var evolution = TestEvolution.init(
        allocator,
        state[0..],
        stepper,
        aux,
    );
    defer evolution.deinit();

    try testing.expectEqual(state.len, evolution.exactValues().len);
    try testing.expectEqual(state.len, evolution.stepper.rhs.len);
    try testing.expectEqual(state.len, evolution.stepper.solution.len);
}

test "Evolution computes error against the current exact field" {
    const TestEvolution = Evolution([]f64, MockStepper, MockExactAux);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };
    const builder = MockBuilder{ .scratch_len = state.len };
    const stepper = try builder.initStepper(allocator, state[0..]);
    const aux = try MockExactAux.init(
        allocator,
        &mesh,
        state.len,
        MockExactInitializer{},
        MockErrorMeasure{},
    );

    var evolution = TestEvolution.init(
        allocator,
        state[0..],
        stepper,
        aux,
    );
    defer evolution.deinit();

    evolution.fillExact(0.0);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(2.0, evolution.l2Error(), 1e-15);
}
