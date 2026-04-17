//! Library-side evolution mechanics for reusable system-stepping workflows.
//!
//! `Evolution` is the execution noun: it owns time-step policy, run counters,
//! and listeners. The owned PDE runtime lives on the chosen system object.
//! Analytic/reference comparisons are separate higher-order abstractions and do
//! not live inside `Evolution`.

const std = @import("std");
const testing = std.testing;

pub const reference = @import("reference_study.zig");

pub fn Method(
    comptime SystemType: type,
    comptime MethodType: type,
) void {
    const MethodOptionsType = methodOptionsType(MethodType);

    if (!@hasDecl(MethodType, "advance")) {
        @compileError("Method requires `pub fn advance(allocator, system, dt[, options]) !void`");
    }

    const advance_info = @typeInfo(@TypeOf(MethodType.advance));
    if (advance_info != .@"fn") {
        @compileError("Method: `advance` must be a function");
    }

    const advance_fn_info = advance_info.@"fn";
    const expected_param_len: usize = if (MethodOptionsType == void) 3 else 4;
    if (advance_fn_info.params.len != expected_param_len) {
        @compileError("Method: `advance` must take (allocator, system, dt[, options])");
    }
    if (advance_fn_info.params[0].type != std.mem.Allocator) {
        @compileError("Method: `advance` allocator parameter must be std.mem.Allocator");
    }
    if (advance_fn_info.params[1].type != SystemType) {
        @compileError("Method: `advance` system parameter must match Evolution system type");
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
            @compileError("Method: `initialize` must take (allocator, system, dt[, options])");
        }
        if (initialize_fn_info.params[0].type != std.mem.Allocator) {
            @compileError("Method: `initialize` allocator parameter must be std.mem.Allocator");
        }
        if (initialize_fn_info.params[1].type != SystemType) {
            @compileError("Method: `initialize` system parameter must match Evolution system type");
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
    if (!@hasDecl(MethodType, "Options")) return void;
    const Options = MethodType.Options;
    if (@TypeOf(Options) != type) {
        @compileError("Method: `Options` must be a type");
    }
    return Options;
}

fn resolvedMethodOptions(
    comptime MethodType: type,
    maybe_options: if (methodOptionsType(MethodType) == void) void else ?methodOptionsType(MethodType),
) methodOptionsType(MethodType) {
    if (methodOptionsType(MethodType) == void) return {};
    if (maybe_options) |options| return options;
    if (@hasDecl(MethodType, "defaultOptions")) return MethodType.defaultOptions();
    @panic("Evolution.Config.init requires `.methodOptions(...)` for methods without default options");
}

fn initializeMethod(
    comptime MethodType: type,
    allocator: std.mem.Allocator,
    system: anytype,
    time_step: f64,
    options: methodOptionsType(MethodType),
) !void {
    if (!@hasDecl(MethodType, "initialize")) return;

    if (methodOptionsType(MethodType) == void) {
        const result = MethodType.initialize(allocator, system, time_step);
        if (@TypeOf(result) == void) return;
        try result;
        return;
    }

    const result = MethodType.initialize(allocator, system, time_step, options);
    if (@TypeOf(result) == void) return;
    try result;
}

fn advanceMethod(
    comptime MethodType: type,
    allocator: std.mem.Allocator,
    system: anytype,
    time_step: f64,
    options: methodOptionsType(MethodType),
) !void {
    if (methodOptionsType(MethodType) == void) {
        try MethodType.advance(allocator, system, time_step);
        return;
    }
    try MethodType.advance(allocator, system, time_step, options);
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

pub fn Evolution(
    comptime SystemType: type,
    comptime MethodType: type,
) type {
    comptime {
        Method(SystemType, MethodType);
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

            pub fn init(
                self: Config,
                allocator: std.mem.Allocator,
                system: SystemType,
            ) !Self {
                const time_step = self.time_step orelse
                    @panic("Evolution.Config.init requires `.dt(...)` before construction");
                const method_options = resolvedMethodOptions(MethodType, self.method_options);

                const evolution = Self{
                    .allocator = allocator,
                    .system = system,
                    .time_step = time_step,
                    .current_time = 0.0,
                    .step_count = 0,
                    .method_options = method_options,
                };
                try initializeMethod(MethodType, allocator, system, time_step, method_options);
                return evolution;
            }
        };

        allocator: std.mem.Allocator,
        system: SystemType,
        time_step: f64,
        current_time: f64,
        step_count: u32,
        method_options: MethodOptionsType,

        pub fn config() Config {
            return .{};
        }

        pub fn deinit(self: *Self) void {
            _ = self;
        }

        pub fn step(self: *Self, allocator: std.mem.Allocator) !void {
            try advanceMethod(MethodType, allocator, self.system, self.time_step, self.method_options);
            self.step_count += 1;
            self.current_time += self.time_step;
        }

        pub fn currentSystem(self: *const Self) SystemType {
            return self.system;
        }

        pub fn timeStep(self: *const Self) f64 {
            return self.time_step;
        }

        pub fn currentTime(self: *const Self) f64 {
            return self.current_time;
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

const MockMethod = struct {
    pub fn initialize(_: std.mem.Allocator, system_values: []f64, _: f64) void {
        for (system_values, 0..) |*value, idx| {
            value.* = 1.0 + @as(f64, @floatFromInt(idx));
        }
    }

    pub fn advance(_: std.mem.Allocator, system_values: []f64, _: f64) !void {
        for (system_values) |*value| {
            value.* += 1.0;
        }
    }
};

const MockAdaptiveMethod = struct {
    pub const Options = struct {
        gain: f64,
    };

    pub fn advance(_: std.mem.Allocator, system_values: []f64, _: f64, options: Options) !void {
        for (system_values) |*value| {
            value.* += options.gain;
        }
    }
};

test "Evolution config owns dt and initializes the chosen system" {
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 0.0, 0.0 };
    var evolution = try Evolution([]f64, MockMethod).config()
        .dt(0.25)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try testing.expectApproxEqAbs(0.25, evolution.timeStep(), 1e-15);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(1.0, system[0], 1e-15);
    try testing.expectApproxEqAbs(2.0, system[1], 1e-15);
    try testing.expectApproxEqAbs(3.0, system[2], 1e-15);

    try evolution.step(allocator);
    try testing.expectEqual(@as(u32, 1), evolution.stepCount());
    try testing.expectApproxEqAbs(2.0, system[0], 1e-15);
    try testing.expectApproxEqAbs(3.0, system[1], 1e-15);
    try testing.expectApproxEqAbs(4.0, system[2], 1e-15);
    try testing.expectApproxEqAbs(0.25, evolution.currentTime(), 1e-15);
}

test "Evolution method options stay method-local" {
    const allocator = testing.allocator;
    var system = [_]f64{0.0};
    var evolution = try Evolution([]f64, MockAdaptiveMethod).config()
        .dt(0.5)
        .methodOptions(.{ .gain = 3.0 })
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try evolution.step(allocator);
    try testing.expectApproxEqAbs(3.0, system[0], 1e-15);
}

test "Evolution run notifies listeners at begin, each step, and end" {
    const TestEvolution = Evolution([]f64, MockMethod);
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 0.0 };
    var evolution = try TestEvolution.config()
        .dt(0.5)
        .init(allocator, system[0..]);
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
