//! Library-side evolution mechanics for time-stepped systems with exact/reference fields.
//!
//! This layer owns the reusable buffers and step orchestration state that sit
//! between assembled operators and higher-level example run loops. It is
//! intentionally narrower than a universal PDE runtime: the first supported
//! shape is a system that evolves a single field against an exact/reference
//! solution while reusing scratch across steps.

const std = @import("std");
const testing = std.testing;

fn validateBuilderType(comptime BuilderType: type, comptime Value: type) void {
    if (!@hasDecl(BuilderType, "Scratch")) {
        @compileError("evolution builder must declare `pub const Scratch`");
    }
    if (!@hasDecl(BuilderType, "Stepper")) {
        @compileError("evolution builder must declare `pub const Stepper`");
    }
    if (!@hasDecl(BuilderType, "initScratch")) {
        @compileError("evolution builder must declare `pub fn initScratch(self, allocator) !Scratch`");
    }
    if (!@hasDecl(BuilderType, "deinitScratch")) {
        @compileError("evolution builder must declare `pub fn deinitScratch(self, allocator, *Scratch) void`");
    }
    if (!@hasDecl(BuilderType, "seedSolution")) {
        @compileError("evolution builder must declare `pub fn seedSolution(self, state_values, *Scratch) void`");
    }
    if (!@hasDecl(BuilderType, "makeStepper")) {
        @compileError("evolution builder must declare `pub fn makeStepper(self, state_values, *Scratch) Stepper`");
    }

    const StepperType = BuilderType.Stepper;
    if (!@hasDecl(StepperType, "step")) {
        @compileError("evolution builder stepper must declare `pub fn step(self, allocator) !void`");
    }

    const step_info = @typeInfo(@TypeOf(StepperType.step));
    if (step_info != .@"fn") {
        @compileError("evolution stepper `step` must be a function");
    }

    _ = Value;
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

/// Evolution object for a single time-stepped field with an exact/reference solution.
///
/// The caller owns the mesh and the primary state storage. The evolution object
/// owns the reusable exact/reference buffer and any scratch requested by the
/// builder. This keeps example code focused on mesh/system construction while
/// the library owns repeated timestep orchestration state.
pub fn ExactEvolution(
    comptime MeshType: type,
    comptime Value: type,
    comptime BuilderType: type,
    comptime ExactInitializerType: type,
    comptime ErrorMeasureType: type,
) type {
    comptime {
        validateBuilderType(BuilderType, Value);
        validateExactInitializerType(MeshType, Value, ExactInitializerType);
        validateErrorMeasureType(MeshType, Value, ErrorMeasureType);
    }

    return struct {
        const Self = @This();
        pub const Stepper = BuilderType.Stepper;
        pub const Scratch = BuilderType.Scratch;

        allocator: std.mem.Allocator,
        mesh: *const MeshType,
        state_values: []Value,
        exact_values: []Value,
        builder: BuilderType,
        scratch: Scratch,
        stepper: Stepper,
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

            var scratch = try builder.initScratch(allocator);
            errdefer builder.deinitScratch(allocator, &scratch);

            builder.seedSolution(state_values, &scratch);

            return .{
                .allocator = allocator,
                .mesh = mesh,
                .state_values = state_values,
                .exact_values = exact_values,
                .builder = builder,
                .scratch = scratch,
                .stepper = builder.makeStepper(state_values, &scratch),
                .exact_initializer = exact_initializer,
                .error_measure = error_measure,
            };
        }

        pub fn deinit(self: *Self) void {
            self.builder.deinitScratch(self.allocator, &self.scratch);
            self.allocator.free(self.exact_values);
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try self.stepper.step(allocator);
        }

        pub fn fillExact(self: *Self, time: f64) void {
            self.exact_initializer.fill(self.mesh, self.exact_values, time);
        }

        pub fn l2Error(self: *const Self) f64 {
            return self.error_measure.compute(self.mesh, self.state_values, self.exact_values);
        }

        pub fn stateValues(self: *const Self) []const Value {
            return self.state_values;
        }

        pub fn exactValues(self: *Self) []Value {
            return self.exact_values;
        }
    };
}

const MockMesh = struct {};

const MockScratch = struct {
    rhs: []f64,
    solution: []f64,
};

const MockStepper = struct {
    state_values: []f64,
    rhs: []f64,
    solution: []f64,

    pub fn step(self: @This(), _: std.mem.Allocator) !void {
        for (self.state_values, self.rhs, self.solution) |*state_value, *rhs_value, *solution_value| {
            rhs_value.* = state_value.*;
            solution_value.* = rhs_value.* + 1.0;
            state_value.* = solution_value.*;
        }
    }
};

const MockBuilder = struct {
    pub const Scratch = MockScratch;
    pub const Stepper = MockStepper;

    scratch_len: usize,

    pub fn initScratch(self: @This(), allocator: std.mem.Allocator) !Scratch {
        return .{
            .rhs = try allocator.alloc(f64, self.scratch_len),
            .solution = try allocator.alloc(f64, self.scratch_len),
        };
    }

    pub fn deinitScratch(_: @This(), allocator: std.mem.Allocator, scratch: *Scratch) void {
        allocator.free(scratch.solution);
        allocator.free(scratch.rhs);
    }

    pub fn seedSolution(_: @This(), state_values: []const f64, scratch: *Scratch) void {
        @memcpy(scratch.solution, state_values);
    }

    pub fn makeStepper(_: @This(), state_values: []f64, scratch: *Scratch) Stepper {
        return .{
            .state_values = state_values,
            .rhs = scratch.rhs,
            .solution = scratch.solution,
        };
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

test "ExactEvolution owns exact buffer and reusable scratch" {
    const Evolution = ExactEvolution(MockMesh, f64, MockBuilder, MockExactInitializer, MockErrorMeasure);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0, 3.0 };

    var evolution = try Evolution.init(
        allocator,
        &mesh,
        &state,
        .{ .scratch_len = state.len },
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer evolution.deinit();

    try testing.expectEqual(state.len, evolution.exactValues().len);
    try testing.expectEqual(state.len, evolution.scratch.rhs.len);
    try testing.expectEqual(state.len, evolution.scratch.solution.len);

    evolution.fillExact(0.5);
    try testing.expectApproxEqAbs(0.5, evolution.exactValues()[0], 1e-15);

    try evolution.step(allocator);
    try testing.expectApproxEqAbs(2.0, state[0], 1e-15);
    try testing.expectApproxEqAbs(3.0, state[1], 1e-15);
    try testing.expectApproxEqAbs(4.0, state[2], 1e-15);
}

test "ExactEvolution computes error against the current exact field" {
    const Evolution = ExactEvolution(MockMesh, f64, MockBuilder, MockExactInitializer, MockErrorMeasure);
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };

    var evolution = try Evolution.init(
        allocator,
        &mesh,
        &state,
        .{ .scratch_len = state.len },
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer evolution.deinit();

    evolution.fillExact(0.0);
    try testing.expectApproxEqAbs(2.0, evolution.l2Error(), 1e-15);
}
