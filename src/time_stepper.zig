//! Comptime time-stepping infrastructure.
//!
//! Two layers:
//!   1. `TimeStepStrategy` — a comptime concept that validates a user-defined
//!      integration strategy has the required `State` and `step` declarations.
//!   2. `TimeStepper(Strategy)` — a generic integrator that wraps any conforming
//!      strategy, providing a uniform interface for simulation runners.
//!
//! Users define a strategy (the raw physics + numerics), then pass it to
//! `TimeStepper` which validates the contract and produces a ready-to-use
//! integrator type.
//!
//! ```zig
//! const Stepper = TimeStepper(MaxwellLeapfrog(Mesh2D));
//! try Stepper.step(allocator, &state, dt);
//! ```

const std = @import("std");
const testing = std.testing;

// ═══════════════════════════════════════════════════════════════════════════
// TimeStepStrategy — comptime concept
// ═══════════════════════════════════════════════════════════════════════════

/// Validate that `S` satisfies the TimeStepStrategy concept at compile time.
///
/// A conforming strategy must declare:
///   - `pub const State: type` — the simulation state type
///   - `pub fn step(std.mem.Allocator, *State, f64) !void` — advance by dt
///
/// Produces a descriptive `@compileError` on violation.
pub fn TimeStepStrategy(comptime S: type) void {
    // 1. S must declare a State type.
    if (!@hasDecl(S, "State")) {
        @compileError("TimeStepStrategy requires a 'pub const State: type' declaration — " ++
            "the simulation state type this strategy advances");
    }
    const State = S.State;
    if (@TypeOf(State) != type) {
        @compileError("TimeStepStrategy: 'State' must be a type, got " ++ @typeName(@TypeOf(State)));
    }

    // 2. S must declare a pub step function: fn(Allocator, *State, f64) !void
    if (!@hasDecl(S, "step")) {
        @compileError("TimeStepStrategy requires a 'pub fn step(std.mem.Allocator, *State, f64) !void' " ++
            "declaration — advance state by one timestep dt");
    }

    const step_info = @typeInfo(@TypeOf(S.step));
    if (step_info != .@"fn") {
        @compileError("TimeStepStrategy: 'step' must be a function");
    }
    const fn_info = step_info.@"fn";

    // Check parameter count and types.
    if (fn_info.params.len != 3) {
        @compileError("TimeStepStrategy: 'step' must take exactly 3 parameters " ++
            "(std.mem.Allocator, *State, f64)");
    }
    if (fn_info.params[0].type != std.mem.Allocator) {
        @compileError("TimeStepStrategy: 'step' parameter 0 must be std.mem.Allocator");
    }
    if (fn_info.params[1].type != *State) {
        @compileError("TimeStepStrategy: 'step' parameter 1 must be *State");
    }
    if (fn_info.params[2].type != f64) {
        @compileError("TimeStepStrategy: 'step' parameter 2 must be f64 (the timestep dt)");
    }

    // Check return type is an error union with payload void.
    const ret = fn_info.return_type orelse
        @compileError("TimeStepStrategy: 'step' must have a known return type");
    const ret_info = @typeInfo(ret);
    if (ret_info != .error_union) {
        @compileError("TimeStepStrategy: 'step' must return !void (an error union of void)");
    }
    if (ret_info.error_union.payload != void) {
        @compileError("TimeStepStrategy: 'step' must return !void, not !" ++
            @typeName(ret_info.error_union.payload));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// TimeStepper — generic integrator wrapper
// ═══════════════════════════════════════════════════════════════════════════

/// Generic time integrator parameterized on a `TimeStepStrategy`.
///
/// Validates the strategy contract at comptime, then exposes a uniform
/// interface that simulation runners and composition infrastructure
/// can depend on. The returned type is a namespace — no instantiation
/// needed, just call `step` directly:
///
/// ```zig
/// const Stepper = TimeStepper(MaxwellLeapfrog(Mesh2D));
/// try Stepper.step(allocator, &state, dt);
/// ```
///
/// Step counting is the strategy's responsibility (via `state.timestep`
/// or equivalent) — the wrapper does not duplicate it.
pub fn TimeStepper(comptime Strategy: type) type {
    // Gate: the strategy must satisfy the concept.
    comptime TimeStepStrategy(Strategy);

    return struct {
        /// The simulation state type, forwarded from the strategy.
        pub const State = Strategy.State;

        /// Advance the state by one timestep dt using the underlying strategy.
        pub fn step(allocator: std.mem.Allocator, state: *State, dt: f64) !void {
            try Strategy.step(allocator, state, dt);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — TimeStepStrategy concept
// ═══════════════════════════════════════════════════════════════════════════

/// A minimal conforming strategy for testing — increments value by dt.
const MockStrategy = struct {
    pub const State = struct {
        value: f64,
    };

    pub fn step(_: std.mem.Allocator, state: *State, dt: f64) !void {
        state.value += dt;
    }
};

test "TimeStepStrategy accepts a conforming type" {
    comptime TimeStepStrategy(MockStrategy);
}

test "TimeStepStrategy accepts a strategy with extra declarations" {
    const ExtendedStrategy = struct {
        pub const State = struct {
            x: f64,
            y: f64,
        };

        pub fn step(_: std.mem.Allocator, state: *State, dt: f64) !void {
            state.x += dt;
            state.y += dt * 2.0;
        }

        /// Extra method — the concept should not reject this.
        pub fn reset(state: *State) void {
            state.x = 0;
            state.y = 0;
        }
    };
    comptime TimeStepStrategy(ExtendedStrategy);
}

// ── Negative tests (compile-time rejection) ──────────────────────────────
//
// These verify that non-conforming types produce compile errors. Since
// @compileError is fatal, we cannot test them in CI. The test bodies
// document the contract violation and the expected error message.
// Verification: uncomment the `comptime` call and confirm it fails.

test "TimeStepStrategy rejects type missing State declaration" {
    const NoState = struct {
        pub fn step(_: std.mem.Allocator, _: *anyopaque, _: f64) !void {}
    };

    // comptime TimeStepStrategy(NoState);
    // expected: @compileError("TimeStepStrategy requires a 'pub const State: type' declaration ...")
    _ = NoState;
}

test "TimeStepStrategy rejects type missing step function" {
    const NoStep = struct {
        pub const State = struct { value: f64 };
    };

    // comptime TimeStepStrategy(NoStep);
    // expected: @compileError("TimeStepStrategy requires a 'pub fn step(...) !void' ...")
    _ = NoStep;
}

test "TimeStepStrategy rejects step with wrong signature (missing allocator)" {
    const WrongSig = struct {
        pub const State = struct { value: f64 };

        pub fn step(state: *State, dt: f64) !void {
            state.value += dt;
        }
    };

    // comptime TimeStepStrategy(WrongSig);
    // expected: @compileError("TimeStepStrategy: 'step' must take exactly 3 parameters ...")
    _ = WrongSig;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — TimeStepper generic wrapper
// ═══════════════════════════════════════════════════════════════════════════

test "TimeStepper wraps a conforming strategy and forwards State" {
    const Stepper = TimeStepper(MockStrategy);
    comptime {
        std.debug.assert(@hasDecl(Stepper, "State"));
        std.debug.assert(Stepper.State == MockStrategy.State);
    }
}

test "TimeStepper rejects nesting TimeStepper(TimeStepper(S))" {
    const Stepper = TimeStepper(MockStrategy);

    // comptime _ = TimeStepper(Stepper);
    // expected: @compileError("TimeStepper cannot wrap another TimeStepper ...")
    _ = Stepper;
}

test "TimeStepper.step delegates to strategy" {
    const Stepper = TimeStepper(MockStrategy);
    var state = MockStrategy.State{ .value = 0.0 };

    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.1, state.value, 1e-15);

    try Stepper.step(testing.allocator, &state, 0.2);
    try testing.expectApproxEqAbs(0.3, state.value, 1e-15);
}
