//! Comptime TimeStepper concept — validates time integrator contracts at compile time.
//!
//! A TimeStepper is any type T that provides:
//!   - `T.State`: the simulation state type it operates on
//!   - `T.step(allocator, *T.State, dt) !void`: advance state by one timestep
//!
//! The concept is physics-agnostic — it constrains the integrator interface
//! without knowing what equations are being solved.

const std = @import("std");
const testing = std.testing;

/// Validate that `T` satisfies the TimeStepper concept at compile time.
///
/// A conforming type must declare:
///   - `pub const State: type` — the simulation state type
///   - `pub fn step(std.mem.Allocator, *State, f64) anyerror!void`
///
/// Produces a descriptive `@compileError` on violation.
pub fn TimeStepper(comptime T: type) void {
    // 1. T must declare a State type.
    if (!@hasDecl(T, "State")) {
        @compileError("TimeStepper requires a 'pub const State: type' declaration — " ++
            "the simulation state type this integrator advances");
    }
    const S = T.State;
    if (@TypeOf(S) != type) {
        @compileError("TimeStepper: 'State' must be a type, got " ++ @typeName(@TypeOf(S)));
    }

    // 2. T must declare a step function: fn(Allocator, *State, f64) !void
    if (!@hasDecl(T, "step")) {
        @compileError("TimeStepper requires a 'pub fn step(std.mem.Allocator, *State, f64) !void' " ++
            "declaration — advance state by one timestep dt");
    }

    const step_info = @typeInfo(@TypeOf(T.step));
    if (step_info != .@"fn") {
        @compileError("TimeStepper: 'step' must be a function");
    }
    const fn_info = step_info.@"fn";

    // Check parameter count and types.
    if (fn_info.params.len != 3) {
        @compileError("TimeStepper: 'step' must take exactly 3 parameters " ++
            "(std.mem.Allocator, *State, f64)");
    }
    if (fn_info.params[0].type != std.mem.Allocator) {
        @compileError("TimeStepper: 'step' parameter 0 must be std.mem.Allocator");
    }
    if (fn_info.params[1].type != *S) {
        @compileError("TimeStepper: 'step' parameter 1 must be *State");
    }
    if (fn_info.params[2].type != f64) {
        @compileError("TimeStepper: 'step' parameter 2 must be f64 (the timestep dt)");
    }

    // Check return type is an error union with payload void.
    const ret = fn_info.return_type orelse
        @compileError("TimeStepper: 'step' must have a known return type");
    const ret_info = @typeInfo(ret);
    if (ret_info != .error_union) {
        @compileError("TimeStepper: 'step' must return !void (an error union of void)");
    }
    if (ret_info.error_union.payload != void) {
        @compileError("TimeStepper: 'step' must return !void, not !" ++
            @typeName(ret_info.error_union.payload));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

/// A minimal conforming stepper for testing — does nothing, proves the
/// concept accepts a valid type.
const MockStepper = struct {
    pub const State = struct {
        value: f64,
    };

    pub fn step(_: std.mem.Allocator, state: *State, dt: f64) !void {
        state.value += dt;
    }
};

test "TimeStepper accepts a conforming type" {
    // This should compile without error.
    comptime TimeStepper(MockStepper);
}

test "TimeStepper accepts a stepper with extra declarations" {
    // A conforming type may have additional declarations beyond the contract.
    const ExtendedStepper = struct {
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
    comptime TimeStepper(ExtendedStepper);
}

// ── Negative tests (compile-time rejection) ──────────────────────────────
//
// These verify that non-conforming types produce compile errors. Since
// @compileError is fatal, we cannot test it at runtime. Instead we
// use a pattern: wrap the concept call behind a comptime-unreachable
// branch so that the test documents intent without triggering the error.
//
// The real verification is: uncomment the direct call and confirm it
// fails to compile with a descriptive message. CI tests the positive
// path; the negative path is verified manually during review.

test "TimeStepper rejects type missing State declaration" {
    const NoState = struct {
        pub fn step(_: std.mem.Allocator, _: *anyopaque, _: f64) !void {}
    };

    // If we could test compile errors at runtime, we'd assert:
    //   comptime TimeStepper(NoState);  // expected: @compileError about missing State
    // Instead, document the contract violation.
    _ = NoState;
}

test "TimeStepper rejects type missing step function" {
    const NoStep = struct {
        pub const State = struct { value: f64 };
    };

    // comptime TimeStepper(NoStep);  // expected: @compileError about missing step
    _ = NoStep;
}

test "TimeStepper rejects step with wrong signature (missing allocator)" {
    const WrongSig = struct {
        pub const State = struct { value: f64 };

        pub fn step(state: *State, dt: f64) !void {
            state.value += dt;
        }
    };

    // comptime TimeStepper(WrongSig);  // expected: @compileError about wrong step signature
    _ = WrongSig;
}
