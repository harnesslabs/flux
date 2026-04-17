//! Library-side evolution mechanics for reusable state-stepping systems.
//!
//! `Evolution` is the execution noun: it owns time-step policy, run counters,
//! listeners, and optional exact/reference auxiliary state. The integration
//! method is a comptime type parameter. Method-specific setup flows through
//! `Evolution.Config`, not through a separate runtime "stepper" object that
//! redundantly owns state or dt.

const std = @import("std");
const testing = std.testing;

pub fn Method(
    comptime StateType: type,
    comptime MethodType: type,
) void {
    const MethodOptionsType = methodOptionsType(MethodType);

    if (!@hasDecl(MethodType, "advance")) {
        @compileError("Method requires `pub fn advance(allocator, state, dt[, options]) !void`");
    }

    const advance_info = @typeInfo(@TypeOf(MethodType.advance));
    if (advance_info != .@"fn") {
        @compileError("Method: `advance` must be a function");
    }

    const advance_fn_info = advance_info.@"fn";
    const expected_param_len: usize = if (MethodOptionsType == void) 3 else 4;
    if (advance_fn_info.params.len != expected_param_len) {
        @compileError("Method: `advance` must take (allocator, state, dt[, options])");
    }
    if (advance_fn_info.params[0].type != std.mem.Allocator) {
        @compileError("Method: `advance` allocator parameter must be std.mem.Allocator");
    }
    if (advance_fn_info.params[1].type != StateType) {
        @compileError("Method: `advance` state parameter must match Evolution state type");
    }
    if (advance_fn_info.params[2].type != f64) {
        @compileError("Method: `advance` timestep parameter must be f64");
    }
    if (MethodOptionsType != void and advance_fn_info.params[3].type != MethodOptionsType) {
        @compileError("Method: `advance` options parameter must match `pub const Options`");
    }

    const advance_return = advance_fn_info.return_type orelse
        @compileError("Method: `advance` must have a known return type");
    const advance_return_info = @typeInfo(advance_return);
    if (advance_return_info != .error_union or advance_return_info.error_union.payload != void) {
        @compileError("Method: `advance` must return !void");
    }

    if (@hasDecl(MethodType, "initialize")) {
        const initialize_info = @typeInfo(@TypeOf(MethodType.initialize));
        if (initialize_info != .@"fn") {
            @compileError("Method: `initialize` must be a function");
        }

        const initialize_fn_info = initialize_info.@"fn";
        const expected_initialize_params: usize = if (MethodOptionsType == void) 3 else 4;
        if (initialize_fn_info.params.len != expected_initialize_params) {
            @compileError("Method: `initialize` must take (allocator, state, dt[, options])");
        }
        if (initialize_fn_info.params[0].type != std.mem.Allocator) {
            @compileError("Method: `initialize` allocator parameter must be std.mem.Allocator");
        }
        if (initialize_fn_info.params[1].type != StateType) {
            @compileError("Method: `initialize` state parameter must match Evolution state type");
        }
        if (initialize_fn_info.params[2].type != f64) {
            @compileError("Method: `initialize` timestep parameter must be f64");
        }
        if (MethodOptionsType != void and initialize_fn_info.params[3].type != MethodOptionsType) {
            @compileError("Method: `initialize` options parameter must match `pub const Options`");
        }

        const initialize_return = initialize_fn_info.return_type orelse
            @compileError("Method: `initialize` must have a known return type");
        if (initialize_return != void) {
            const initialize_return_info = @typeInfo(initialize_return);
            if (initialize_return_info != .error_union or initialize_return_info.error_union.payload != void) {
                @compileError("Method: `initialize` must return void or !void");
            }
        }
    }
}

fn methodOptionsType(comptime MethodType: type) type {
    if (@hasDecl(MethodType, "Options")) {
        const Options = MethodType.Options;
        if (@TypeOf(Options) != type) {
            @compileError("Method: `Options` must be a type");
        }
        return Options;
    }
    return void;
}

fn defaultMethodOptions(comptime MethodType: type) methodOptionsType(MethodType) {
    const MethodOptionsType = methodOptionsType(MethodType);
    if (MethodOptionsType == void) return {};
    if (@hasDecl(MethodType, "defaultOptions")) {
        return MethodType.defaultOptions();
    }
    return .{};
}

fn resolvedMethodOptions(
    comptime MethodType: type,
    maybe_options: if (methodOptionsType(MethodType) == void) void else ?methodOptionsType(MethodType),
) methodOptionsType(MethodType) {
    if (methodOptionsType(MethodType) == void) return {};
    if (maybe_options) |options| return options;
    if (@hasDecl(MethodType, "defaultOptions")) {
        return MethodType.defaultOptions();
    }
    @panic("Evolution.Config.init requires `.methodOptions(...)` for methods without default options");
}

fn requireOptionField(
    comptime MethodType: type,
    comptime field_name: []const u8,
) void {
    const MethodOptionsType = methodOptionsType(MethodType);
    if (MethodOptionsType == void or !@hasField(MethodOptionsType, field_name)) {
        @compileError("Evolution.Config: chosen method does not support `" ++ field_name ++ "`");
    }
}

fn setOptionField(
    comptime MethodType: type,
    options: methodOptionsType(MethodType),
    comptime field_name: []const u8,
    value: anytype,
) methodOptionsType(MethodType) {
    requireOptionField(MethodType, field_name);
    var next = options;
    @field(next, field_name) = value;
    return next;
}

fn invokeListeners(
    comptime method_name: []const u8,
    listeners: anytype,
    event: anytype,
) !void {
    inline for (listeners) |listener| {
        const ListenerDeclType = switch (@typeInfo(@TypeOf(listener))) {
            .pointer => |pointer| pointer.child,
            else => @TypeOf(listener),
        };
        if (@hasDecl(ListenerDeclType, method_name)) {
            try @field(ListenerDeclType, method_name)(listener, event);
        }
    }
}

fn maybeInitializeMethod(
    comptime StateType: type,
    comptime MethodType: type,
    allocator: std.mem.Allocator,
    state: StateType,
    time_step: f64,
    options: methodOptionsType(MethodType),
) !void {
    _ = comptime StateType;
    if (!@hasDecl(MethodType, "initialize")) return;

    if (methodOptionsType(MethodType) == void) {
        const result = MethodType.initialize(allocator, state, time_step);
        const ResultType = @TypeOf(result);
        if (ResultType == void) return;
        try result;
        return;
    }

    const result = MethodType.initialize(allocator, state, time_step, options);
    const ResultType = @TypeOf(result);
    if (ResultType == void) return;
    try result;
}

fn advanceMethod(
    comptime StateType: type,
    comptime MethodType: type,
    allocator: std.mem.Allocator,
    state: StateType,
    time_step: f64,
    options: methodOptionsType(MethodType),
) !void {
    _ = comptime StateType;
    if (methodOptionsType(MethodType) == void) {
        try MethodType.advance(allocator, state, time_step);
        return;
    }
    try MethodType.advance(allocator, state, time_step, options);
}

pub fn ReferenceAux(
    comptime MeshType: type,
    comptime InitializerType: type,
    comptime ErrorMeasureType: type,
) type {
    return struct {
        mesh: *const MeshType,
        exact_values: []f64,
        initializer: InitializerType,
        error_measure: ErrorMeasureType,

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            len: usize,
            initializer: InitializerType,
            error_measure: ErrorMeasureType,
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
}

pub fn empiricalRates(
    allocator: std.mem.Allocator,
    errors: []const f64,
    refinement_ratio: f64,
) ![]f64 {
    std.debug.assert(refinement_ratio > 1.0);

    if (errors.len < 2) {
        return allocator.alloc(f64, 0);
    }

    const rates = try allocator.alloc(f64, errors.len - 1);
    errdefer allocator.free(rates);

    for (0..errors.len - 1) |idx| {
        const coarse = errors[idx];
        const fine = errors[idx + 1];
        std.debug.assert(coarse > 0.0);
        std.debug.assert(fine > 0.0);
        rates[idx] = std.math.log(f64, refinement_ratio, coarse / fine);
    }

    return rates;
}

/// Evolution object for reusable state stepping.
///
/// The caller owns the primary state storage. `Evolution` owns execution state:
/// timestep, current time, step count, listeners, and optional exact/reference
/// auxiliary state. The method type is a comptime parameter and method-specific
/// configuration flows through `Config`.
pub fn Evolution(
    comptime StateType: type,
    comptime MethodType: type,
    comptime AuxType: type,
) type {
    comptime {
        Method(StateType, MethodType);
    }

    const MethodOptionsType = methodOptionsType(MethodType);

    return struct {
        const Self = @This();

        pub const RunConfig = struct {
            steps: u32,
            final_time: f64,
        };

        pub const RunResult = struct {
            elapsed_s: f64,
        };

        pub const Event = struct {
            evolution: *Self,
            step_index: u32,
            time: f64,
            final_time: f64,
        };

        pub const Config = struct {
            time_step: ?f64 = null,
            method_options: if (MethodOptionsType == void) void else ?MethodOptionsType = if (MethodOptionsType == void) {} else null,

            pub fn dt(self: Config, value: f64) Config {
                std.debug.assert(value > 0.0);
                var next = self;
                next.time_step = value;
                return next;
            }

            pub fn methodOptions(self: Config, value: MethodOptionsType) Config {
                if (MethodOptionsType == void) {
                    @compileError("Evolution.Config: chosen method does not expose `Options`");
                }
                var next = self;
                next.method_options = value;
                return next;
            }

            pub fn minStepSize(self: Config, value: f64) Config {
                std.debug.assert(value > 0.0);
                var next = self;
                next.method_options = setOptionField(
                    MethodType,
                    resolvedMethodOptions(MethodType, self.method_options),
                    "min_step_size",
                    value,
                );
                return next;
            }

            pub fn maxStepSize(self: Config, value: f64) Config {
                std.debug.assert(value > 0.0);
                var next = self;
                next.method_options = setOptionField(
                    MethodType,
                    resolvedMethodOptions(MethodType, self.method_options),
                    "max_step_size",
                    value,
                );
                return next;
            }

            pub fn tolerance(self: Config, value: f64) Config {
                std.debug.assert(value > 0.0);
                var next = self;
                next.method_options = setOptionField(
                    MethodType,
                    resolvedMethodOptions(MethodType, self.method_options),
                    "tolerance",
                    value,
                );
                return next;
            }

            pub fn init(
                self: Config,
                allocator: std.mem.Allocator,
                state: StateType,
                aux: AuxType,
            ) !Self {
                const time_step = self.time_step orelse
                    @panic("Evolution.Config.init requires `.dt(...)` before construction");
                const method_options = resolvedMethodOptions(MethodType, self.method_options);

                const evolution = Self{
                    .allocator = allocator,
                    .state = state,
                    .aux = aux,
                    .time_step = time_step,
                    .current_time = 0.0,
                    .step_count = 0,
                    .method_options = method_options,
                };
                try maybeInitializeMethod(StateType, MethodType, allocator, state, time_step, method_options);
                return evolution;
            }
        };

        allocator: std.mem.Allocator,
        state: StateType,
        aux: AuxType,
        time_step: f64,
        current_time: f64,
        step_count: u32,
        method_options: MethodOptionsType,

        pub fn config() Config {
            return .{};
        }

        pub fn deinit(self: *Self) void {
            if (AuxType != void) {
                self.aux.deinit(self.allocator);
            }
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try advanceMethod(StateType, MethodType, allocator, self.state, self.time_step, self.method_options);
            self.step_count += 1;
            self.current_time += self.time_step;
        }

        pub fn currentState(self: *const Self) StateType {
            return self.state;
        }

        pub fn stateValues(self: *const Self) StateType {
            return self.state;
        }

        pub fn timeStep(self: *const Self) f64 {
            return self.time_step;
        }

        pub fn currentTime(self: *const Self) f64 {
            return self.current_time;
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

        pub fn run(self: *Self, allocator: std.mem.Allocator, run_config: RunConfig, listeners: anytype) !RunResult {
            std.debug.assert(run_config.steps > 0);
            std.debug.assert(self.time_step > 0.0);

            const expected_final_time = self.current_time + self.time_step * @as(f64, @floatFromInt(run_config.steps));
            std.debug.assert(std.math.approxEqAbs(f64, expected_final_time, run_config.final_time, 1e-12) or
                std.math.approxEqRel(f64, expected_final_time, run_config.final_time, 1e-12));

            try invokeListeners("onRunBegin", listeners, Event{
                .evolution = self,
                .step_index = self.step_count,
                .time = self.current_time,
                .final_time = run_config.final_time,
            });

            const start_ns = std.time.nanoTimestamp();
            for (0..run_config.steps) |_| {
                try self.step(allocator);
                try invokeListeners("onStep", listeners, Event{
                    .evolution = self,
                    .step_index = self.step_count,
                    .time = self.current_time,
                    .final_time = run_config.final_time,
                });
            }
            const elapsed_ns = std.time.nanoTimestamp() - start_ns;

            try invokeListeners("onRunEnd", listeners, Event{
                .evolution = self,
                .step_index = self.step_count,
                .time = self.current_time,
                .final_time = run_config.final_time,
            });

            return .{
                .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
            };
        }
    };
}

const MockMesh = struct {};

const MockMethod = struct {
    pub const Options = struct {
        increment: f64 = 1.0,
    };

    pub fn defaultOptions() Options {
        return .{};
    }

    pub fn initialize(_: std.mem.Allocator, state_values: []f64, _: f64, _: Options) void {
        for (state_values, 0..) |*value, idx| {
            value.* = 1.0 + @as(f64, @floatFromInt(idx));
        }
    }

    pub fn advance(_: std.mem.Allocator, state_values: []f64, _: f64, options: Options) !void {
        for (state_values) |*state_value| {
            state_value.* += options.increment;
        }
    }
};

const MockAdaptiveMethod = struct {
    pub const Options = struct {
        min_step_size: f64 = 1.0e-4,
        max_step_size: f64 = 1.0,
        tolerance: f64 = 1.0e-6,
    };

    pub fn defaultOptions() Options {
        return .{};
    }

    pub fn advance(_: std.mem.Allocator, state_values: []f64, dt: f64, options: Options) !void {
        const bounded_dt = @min(options.max_step_size, @max(options.min_step_size, dt));
        for (state_values) |*state_value| {
            state_value.* += bounded_dt + options.tolerance;
        }
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

test "Evolution config owns dt and method options" {
    const allocator = testing.allocator;
    var state = [_]f64{ 0.0, 0.0, 0.0 };
    var evolution = try Evolution([]f64, MockMethod, void).config()
        .dt(0.25)
        .methodOptions(.{ .increment = 2.0 })
        .init(allocator, state[0..], {});
    defer evolution.deinit();

    try testing.expectApproxEqAbs(0.25, evolution.timeStep(), 1e-15);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(1.0, state[0], 1e-15);
    try testing.expectApproxEqAbs(2.0, state[1], 1e-15);
    try testing.expectApproxEqAbs(3.0, state[2], 1e-15);

    try evolution.step(allocator);
    try testing.expectEqual(@as(u32, 1), evolution.stepCount());
    try testing.expectApproxEqAbs(3.0, state[0], 1e-15);
    try testing.expectApproxEqAbs(4.0, state[1], 1e-15);
    try testing.expectApproxEqAbs(5.0, state[2], 1e-15);
    try testing.expectApproxEqAbs(0.25, evolution.currentTime(), 1e-15);
}

test "Evolution config forwards adaptive method settings through builder methods" {
    const allocator = testing.allocator;
    var state = [_]f64{0.0};
    var evolution = try Evolution([]f64, MockAdaptiveMethod, void).config()
        .dt(2.0)
        .minStepSize(0.2)
        .maxStepSize(0.5)
        .tolerance(1.0e-3)
        .init(allocator, state[0..], {});
    defer evolution.deinit();

    try evolution.step(allocator);
    try testing.expectApproxEqAbs(0.501, state[0], 1e-15);
}

test "Evolution owns exact buffer through auxiliary state" {
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0, 3.0 };
    const aux = try MockExactAux.init(
        allocator,
        &mesh,
        state.len,
        MockExactInitializer{},
        MockErrorMeasure{},
    );

    var evolution = try Evolution([]f64, MockMethod, MockExactAux).config()
        .dt(0.1)
        .init(allocator, state[0..], aux);
    defer evolution.deinit();

    try testing.expectEqual(state.len, evolution.exactValues().len);
}

test "Evolution computes error against the current exact field" {
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };
    const aux = try MockExactAux.init(
        allocator,
        &mesh,
        state.len,
        MockExactInitializer{},
        MockErrorMeasure{},
    );

    var evolution = try Evolution([]f64, MockMethod, MockExactAux).config()
        .dt(0.1)
        .init(allocator, state[0..], aux);
    defer evolution.deinit();

    evolution.fillExact(0.0);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(2.0, evolution.l2Error(), 1e-15);
}

test "ReferenceAux owns exact buffer and computes error generically" {
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };

    const Aux = ReferenceAux(MockMesh, MockExactInitializer, MockErrorMeasure);
    var aux = try Aux.init(
        allocator,
        &mesh,
        state.len,
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer aux.deinit(allocator);

    aux.fillExact(0.0);
    try testing.expectEqual(state.len, aux.exactValues().len);
    try testing.expectApproxEqAbs(2.0, aux.l2Error(state[0..]), 1e-15);
}

test "empiricalRates returns pairwise convergence rates" {
    const allocator = testing.allocator;
    const errors = [_]f64{ 4.0, 1.0, 0.25 };
    const rates = try empiricalRates(allocator, &errors, 2.0);
    defer allocator.free(rates);

    try testing.expectEqual(@as(usize, 2), rates.len);
    try testing.expectApproxEqAbs(2.0, rates[0], 1e-12);
    try testing.expectApproxEqAbs(2.0, rates[1], 1e-12);
}

test "Evolution run notifies listeners at begin, each step, and end" {
    const TestEvolution = Evolution([]f64, MockMethod, void);
    const allocator = testing.allocator;
    var state = [_]f64{ 0.0, 0.0 };
    var evolution = try TestEvolution.config()
        .dt(0.5)
        .init(allocator, state[0..], {});
    defer evolution.deinit();

    const Listener = struct {
        begin_count: u32 = 0,
        step_count: u32 = 0,
        end_count: u32 = 0,
        last_step_index: u32 = 0,
        last_time: f64 = 0.0,

        pub fn onRunBegin(self: *@This(), event: TestEvolution.Event) !void {
            self.begin_count += 1;
            try testing.expectEqual(@as(u32, 0), event.step_index);
            try testing.expectApproxEqAbs(0.0, event.time, 1e-15);
        }

        pub fn onStep(self: *@This(), event: TestEvolution.Event) !void {
            self.step_count += 1;
            self.last_step_index = event.step_index;
            self.last_time = event.time;
        }

        pub fn onRunEnd(self: *@This(), event: TestEvolution.Event) !void {
            self.end_count += 1;
            try testing.expectEqual(@as(u32, 2), event.step_index);
            try testing.expectApproxEqAbs(1.0, event.time, 1e-15);
        }
    };

    var listener = Listener{};
    const result = try evolution.run(allocator, .{
        .steps = 2,
        .final_time = 1.0,
    }, .{&listener});

    try testing.expectEqual(@as(u32, 1), listener.begin_count);
    try testing.expectEqual(@as(u32, 2), listener.step_count);
    try testing.expectEqual(@as(u32, 1), listener.end_count);
    try testing.expectEqual(@as(u32, 2), listener.last_step_index);
    try testing.expectApproxEqAbs(1.0, listener.last_time, 1e-15);
    try testing.expect(result.elapsed_s >= 0.0);
}
