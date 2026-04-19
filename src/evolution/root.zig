//! Library-side evolution mechanics for reusable system-stepping workflows.
//!
//! `Evolution` is the execution noun: it owns time-step policy, run counters,
//! configured listeners, and method-local options. The owned PDE runtime lives
//! on the chosen system object. Analytic/reference comparisons remain separate
//! higher-order helpers layered around execution.

const std = @import("std");
const testing = std.testing;

pub const listeners = @import("listeners.zig");
pub const reference = @import("reference_study.zig");

fn resolveMethodType(comptime SystemType: type, comptime MethodSpec: anytype) type {
    if (@TypeOf(MethodSpec) == type) return MethodSpec;

    const spec_info = @typeInfo(@TypeOf(MethodSpec));
    if (spec_info != .@"fn") {
        @compileError("Evolution requires a method type or a method-family function");
    }

    const Resolved = MethodSpec(SystemType);
    if (@TypeOf(Resolved) != type) {
        @compileError("Evolution method-family must return a type when specialized on the chosen system");
    }
    return Resolved;
}

pub fn Method(
    comptime SystemType: type,
    comptime MethodSpec: anytype,
) void {
    const MethodType = resolveMethodType(SystemType, MethodSpec);
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
    listeners_ptr: anytype,
    event: anytype,
) !void {
    inline for (std.meta.fields(@TypeOf(listeners_ptr.*))) |field_info| {
        const listener_slot = &@field(listeners_ptr.*, field_info.name);
        const SlotType = @TypeOf(listener_slot.*);
        switch (@typeInfo(SlotType)) {
            .pointer => |pointer| {
                const ListenerDeclType = pointer.child;
                if (@hasDecl(ListenerDeclType, method_name)) {
                    try @field(ListenerDeclType, method_name)(listener_slot.*, event);
                }
            },
            else => {
                const ListenerDeclType = SlotType;
                if (@hasDecl(ListenerDeclType, method_name)) {
                    try @field(ListenerDeclType, method_name)(listener_slot, event);
                }
            },
        }
    }
}

pub fn Evolution(
    comptime SystemType: type,
    comptime MethodSpec: anytype,
) type {
    const MethodType = resolveMethodType(SystemType, MethodSpec);
    comptime {
        Method(SystemType, MethodSpec);
    }

    const MethodOptionsType = methodOptionsType(MethodType);

    return struct {
        const Namespace = @This();

        pub const RunResult = struct {
            elapsed_s: f64,
        };

        pub const Event = struct {
            allocator: std.mem.Allocator,
            system: SystemType,
            step_index: u32,
            step_goal: u32,
            time: f64,
            final_time: f64,
        };

        pub fn Config(comptime ListenersType: type) type {
            return struct {
                const Self = @This();

                time_step: ?f64 = null,
                total_steps: ?u32 = null,
                listeners_value: ListenersType,
                method_options: if (MethodOptionsType == void) void else ?MethodOptionsType = if (MethodOptionsType == void) {} else null,

                pub fn dt(self: Self, value: f64) Self {
                    std.debug.assert(value > 0.0);
                    var next = self;
                    next.time_step = value;
                    return next;
                }

                pub fn steps(self: Self, value: u32) Self {
                    std.debug.assert(value > 0);
                    var next = self;
                    next.total_steps = value;
                    return next;
                }

                pub fn listen(self: Self, listener: anytype) Config(@TypeOf(self.listeners_value ++ .{listener})) {
                    return .{
                        .time_step = self.time_step,
                        .total_steps = self.total_steps,
                        .listeners_value = self.listeners_value ++ .{listener},
                        .method_options = self.method_options,
                    };
                }

                pub fn methodOptions(self: Self, value: MethodOptionsType) Self {
                    if (MethodOptionsType == void) {
                        @compileError("Evolution.Config: chosen method does not expose `Options`");
                    }
                    var next = self;
                    next.method_options = value;
                    return next;
                }

                pub fn init(
                    self: Self,
                    allocator: std.mem.Allocator,
                    system: SystemType,
                ) !Runtime(ListenersType) {
                    const time_step = self.time_step orelse
                        @panic("Evolution.Config.init requires `.dt(...)` before construction");
                    const total_steps = self.total_steps orelse
                        @panic("Evolution.Config.init requires `.steps(...)` before construction");
                    const method_options = resolvedMethodOptions(MethodType, self.method_options);

                    const evolution = Runtime(ListenersType){
                        .allocator = allocator,
                        .system = system,
                        .time_step = time_step,
                        .total_steps = total_steps,
                        .current_time = 0.0,
                        .step_count = 0,
                        .listeners_value = self.listeners_value,
                        .method_options = method_options,
                    };
                    try initializeMethod(MethodType, allocator, system, time_step, method_options);
                    return evolution;
                }
            };
        }

        pub fn Runtime(comptime ListenersType: type) type {
            return struct {
                const Self = @This();
                pub const EvolutionNamespace = Namespace;

                allocator: std.mem.Allocator,
                system: SystemType,
                time_step: f64,
                total_steps: u32,
                current_time: f64,
                step_count: u32,
                listeners_value: ListenersType,
                method_options: MethodOptionsType,

                pub fn deinit(self: *Self) void {
                    _ = self;
                }

                pub fn step(self: *Self) !void {
                    std.debug.assert(self.step_count < self.total_steps);
                    try advanceMethod(MethodType, self.allocator, self.system, self.time_step, self.method_options);
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

                pub fn configuredSteps(self: *const Self) u32 {
                    return self.total_steps;
                }

                fn remainingSteps(self: *const Self) u32 {
                    std.debug.assert(self.step_count <= self.total_steps);
                    return self.total_steps - self.step_count;
                }

                fn finalTime(self: *const Self) f64 {
                    return self.current_time + self.time_step * @as(f64, @floatFromInt(self.remainingSteps()));
                }

                fn event(self: *const Self, final_time: f64) Event {
                    return .{
                        .allocator = self.allocator,
                        .system = self.system,
                        .step_index = self.step_count,
                        .step_goal = self.total_steps,
                        .time = self.current_time,
                        .final_time = final_time,
                    };
                }

                fn runInternal(self: *Self, extra_listeners_ptr: anytype) !RunResult {
                    const final_time = self.finalTime();
                    const remaining_steps = self.remainingSteps();
                    var configured_listeners = self.listeners_value;

                    try invokeListeners("onRunBegin", &configured_listeners, self.event(final_time));
                    try invokeListeners("onRunBegin", extra_listeners_ptr, self.event(final_time));

                    const start_ns = std.time.nanoTimestamp();
                    for (0..remaining_steps) |_| {
                        try self.step();
                        try invokeListeners("onStep", &configured_listeners, self.event(final_time));
                        try invokeListeners("onStep", extra_listeners_ptr, self.event(final_time));
                    }
                    const elapsed_ns = std.time.nanoTimestamp() - start_ns;

                    try invokeListeners("onRunEnd", &configured_listeners, self.event(final_time));
                    try invokeListeners("onRunEnd", extra_listeners_ptr, self.event(final_time));

                    self.listeners_value = configured_listeners;

                    return .{
                        .elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0,
                    };
                }

                pub fn run(self: *Self) !RunResult {
                    var no_listeners = .{};
                    return self.runInternal(&no_listeners);
                }

                pub fn runWith(self: *Self, extra_listeners: anytype) !RunResult {
                    var listeners_value = extra_listeners;
                    return self.runInternal(&listeners_value);
                }
            };
        }

        pub fn config() Config(@TypeOf(.{})) {
            return .{
                .listeners_value = .{},
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

const MockMethodFamily = struct {
    pub fn specialize(comptime SystemType: type) type {
        return struct {
            pub fn advance(_: std.mem.Allocator, system_values: SystemType, _: f64) !void {
                for (system_values) |*value| {
                    value.* += 2.0;
                }
            }
        };
    }
}.specialize;

test "Evolution config owns dt and initializes the chosen system" {
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 0.0, 0.0 };
    var evolution = try Evolution([]f64, MockMethod).config()
        .dt(0.25)
        .steps(2)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try testing.expectApproxEqAbs(0.25, evolution.timeStep(), 1e-15);
    try testing.expectEqual(@as(u32, 0), evolution.stepCount());
    try testing.expectApproxEqAbs(1.0, system[0], 1e-15);
    try testing.expectApproxEqAbs(2.0, system[1], 1e-15);
    try testing.expectApproxEqAbs(3.0, system[2], 1e-15);

    try evolution.step();
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
        .steps(1)
        .methodOptions(.{ .gain = 3.0 })
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try evolution.step();
    try testing.expectApproxEqAbs(3.0, system[0], 1e-15);
}

test "Evolution specializes method families internally" {
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 1.0 };
    var evolution = try Evolution([]f64, MockMethodFamily).config()
        .dt(1.0)
        .steps(1)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try evolution.step();
    try testing.expectApproxEqAbs(2.0, system[0], 1e-15);
    try testing.expectApproxEqAbs(3.0, system[1], 1e-15);
}

test "Evolution run notifies configured listeners at begin, each step, and end" {
    const EvolutionNs = Evolution([]f64, MockMethod);
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 0.0 };

    const Listener = struct {
        begin_count: u32 = 0,
        step_count: u32 = 0,
        end_count: u32 = 0,
        last_step_index: u32 = 0,
        last_time: f64 = 0.0,

        pub fn onRunBegin(self: *@This(), event: EvolutionNs.Event) !void {
            self.begin_count += 1;
            try testing.expectEqual(@as(u32, 0), event.step_index);
            try testing.expectEqual(@as(u32, 2), event.step_goal);
            try testing.expectApproxEqAbs(0.0, event.time, 1e-15);
        }

        pub fn onStep(self: *@This(), event: EvolutionNs.Event) !void {
            self.step_count += 1;
            self.last_step_index = event.step_index;
            self.last_time = event.time;
        }

        pub fn onRunEnd(self: *@This(), event: EvolutionNs.Event) !void {
            self.end_count += 1;
            try testing.expectEqual(@as(u32, 2), event.step_index);
            try testing.expectApproxEqAbs(1.0, event.time, 1e-15);
        }
    };

    var listener = Listener{};
    var evolution = try EvolutionNs.config()
        .dt(0.5)
        .steps(2)
        .listen(&listener)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    const result = try evolution.run();
    try testing.expect(result.elapsed_s >= 0.0);
    try testing.expectEqual(@as(u32, 1), listener.begin_count);
    try testing.expectEqual(@as(u32, 2), listener.step_count);
    try testing.expectEqual(@as(u32, 1), listener.end_count);
    try testing.expectEqual(@as(u32, 2), listener.last_step_index);
    try testing.expectApproxEqAbs(1.0, listener.last_time, 1e-15);
}

test "Evolution run stops at the configured total after manual stepping" {
    const allocator = testing.allocator;
    var system = [_]f64{0.0};

    var evolution = try Evolution([]f64, MockMethod).config()
        .dt(0.5)
        .steps(2)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    try evolution.step();
    try testing.expectEqual(@as(u32, 1), evolution.stepCount());
    try testing.expectApproxEqAbs(0.5, evolution.currentTime(), 1e-15);

    _ = try evolution.run();

    try testing.expectEqual(@as(u32, 2), evolution.stepCount());
    try testing.expectApproxEqAbs(1.0, evolution.currentTime(), 1e-15);
    try testing.expectApproxEqAbs(3.0, system[0], 1e-15);
}

test "Evolution runWith keeps compatibility for external listeners" {
    const EvolutionNs = Evolution([]f64, MockMethod);
    const allocator = testing.allocator;
    var system = [_]f64{ 0.0, 0.0 };
    var evolution = try EvolutionNs.config()
        .dt(0.5)
        .steps(1)
        .init(allocator, system[0..]);
    defer evolution.deinit();

    const Listener = struct {
        count: u32 = 0,

        pub fn onStep(self: *@This(), _: EvolutionNs.Event) !void {
            self.count += 1;
        }
    };

    var listener = Listener{};
    _ = try evolution.runWith(.{&listener});
    try testing.expectEqual(@as(u32, 1), listener.count);
}
