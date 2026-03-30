//! Generic forward Euler (explicit Euler) integrator.
//!
//! Forward Euler is the simplest explicit time integration scheme:
//!
//!   y^{n+1} = y^n + dt · f(y^n)
//!
//! First-order accurate, unconditionally unstable for oscillatory systems.
//! Useful as a baseline, for prototyping, and for dissipative systems with
//! sufficiently small timesteps.
//!
//! ## Usage
//!
//! Define a system type with a single explicit step, then wrap it:
//!
//! ```zig
//! const MySystem = struct {
//!     pub const State = MyState;
//!     pub fn forward_step(alloc: Allocator, state: *State, dt: f64) !void {
//!         // Compute f(state) and advance: state += dt * f(state)
//!     }
//! };
//! const Stepper = ForwardEuler(MySystem); // satisfies TimeStepStrategy
//! try Stepper.step(allocator, &state, dt);
//! ```

const std = @import("std");
const testing = std.testing;
const time_stepper = @import("../time_stepper.zig");

// ═══════════════════════════════════════════════════════════════════════════
// ExplicitSystem — comptime concept
// ═══════════════════════════════════════════════════════════════════════════

/// Validate that `S` satisfies the ExplicitSystem concept at compile time.
///
/// A conforming system must declare:
///   - `pub const State: type`
///   - `pub fn forward_step(Allocator, *State, f64) !void`
///
/// The `forward_step` computes the RHS and advances the state in one call:
/// `state += dt * f(state)`. This is intentionally simple — the system owns
/// both the derivative evaluation and the state update, giving it full control
/// over intermediate allocations and in-place mutation.
pub fn ExplicitSystem(comptime S: type) void {
    // 1. S must declare a State type.
    if (!@hasDecl(S, "State")) {
        @compileError("ExplicitSystem requires a 'pub const State: type' declaration — " ++
            "the simulation state type this system advances");
    }
    const State = S.State;
    if (@TypeOf(State) != type) {
        @compileError("ExplicitSystem: 'State' must be a type, got " ++ @typeName(@TypeOf(State)));
    }

    // 2. Validate forward_step function signature.
    if (!@hasDecl(S, "forward_step")) {
        @compileError("ExplicitSystem requires a 'pub fn forward_step(std.mem.Allocator, *State, f64) !void' " ++
            "declaration — compute f(state) and advance by dt");
    }

    const info = @typeInfo(@TypeOf(S.forward_step));
    if (info != .@"fn") {
        @compileError("ExplicitSystem: 'forward_step' must be a function");
    }
    const fn_info = info.@"fn";

    if (fn_info.params.len != 3) {
        @compileError("ExplicitSystem: 'forward_step' must take exactly 3 parameters " ++
            "(std.mem.Allocator, *State, f64)");
    }
    if (fn_info.params[0].type != std.mem.Allocator) {
        @compileError("ExplicitSystem: 'forward_step' parameter 0 must be std.mem.Allocator");
    }
    if (fn_info.params[1].type != *State) {
        @compileError("ExplicitSystem: 'forward_step' parameter 1 must be *State");
    }
    if (fn_info.params[2].type != f64) {
        @compileError("ExplicitSystem: 'forward_step' parameter 2 must be f64 (the timestep dt)");
    }

    const ret = fn_info.return_type orelse
        @compileError("ExplicitSystem: 'forward_step' must have a known return type");
    const ret_info = @typeInfo(ret);
    if (ret_info != .error_union) {
        @compileError("ExplicitSystem: 'forward_step' must return !void");
    }
    if (ret_info.error_union.payload != void) {
        @compileError("ExplicitSystem: 'forward_step' must return !void, not !" ++
            @typeName(ret_info.error_union.payload));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ForwardEuler — generic integrator
// ═══════════════════════════════════════════════════════════════════════════

/// Generic forward Euler integrator parameterized on an `ExplicitSystem`.
///
/// Delegates to the system's `forward_step` for a single explicit update.
/// The returned type satisfies `TimeStepStrategy`.
pub fn ForwardEuler(comptime System: type) type {
    comptime ExplicitSystem(System);

    return struct {
        /// The simulation state type, forwarded from the system.
        pub const State = System.State;

        /// Advance the state by one forward Euler step.
        pub fn step(allocator: std.mem.Allocator, state: *State, dt: f64) !void {
            try System.forward_step(allocator, state, dt);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — ExplicitSystem concept
// ═══════════════════════════════════════════════════════════════════════════

/// Mock system: exponential decay dy/dt = -y.
const MockExplicitSystem = struct {
    pub const State = struct {
        value: f64,
    };

    /// Forward step: y += dt * (-y) = y * (1 - dt).
    pub fn forward_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
        state.value += dt * (-state.value);
    }
};

test "ExplicitSystem accepts a conforming type" {
    comptime ExplicitSystem(MockExplicitSystem);
}

test "ExplicitSystem accepts a type with extra declarations" {
    const Extended = struct {
        pub const State = struct { x: f64, y: f64 };
        pub fn forward_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
            const dx = -state.y;
            const dy = state.x;
            state.x += dt * dx;
            state.y += dt * dy;
        }
        /// Extra — concept should not reject.
        pub fn norm(state: *const State) f64 {
            return @sqrt(state.x * state.x + state.y * state.y);
        }
    };
    comptime ExplicitSystem(Extended);
}

// ── Negative tests (compile-time rejection) ──────────────────────────────

test "ExplicitSystem rejects type missing State" {
    const NoState = struct {
        pub fn forward_step(_: std.mem.Allocator, _: *anyopaque, _: f64) !void {}
    };
    // comptime ExplicitSystem(NoState);
    // expected: @compileError("ExplicitSystem requires a 'pub const State: type' declaration ...")
    _ = NoState;
}

test "ExplicitSystem rejects type missing forward_step" {
    const NoStep = struct {
        pub const State = struct { value: f64 };
    };
    // comptime ExplicitSystem(NoStep);
    // expected: @compileError("ExplicitSystem requires a 'pub fn forward_step(...)' ...")
    _ = NoStep;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — ForwardEuler integrator
// ═══════════════════════════════════════════════════════════════════════════

test "ForwardEuler satisfies TimeStepStrategy" {
    const Stepper = ForwardEuler(MockExplicitSystem);
    comptime time_stepper.TimeStepStrategy(Stepper);
}

test "ForwardEuler advances state by one explicit step" {
    const Stepper = ForwardEuler(MockExplicitSystem);
    var state = MockExplicitSystem.State{ .value = 1.0 };

    // dy/dt = -y with y₀ = 1, dt = 0.1:
    //   y₁ = 1.0 + 0.1 * (-1.0) = 0.9
    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.9, state.value, 1e-15);

    // y₂ = 0.9 + 0.1 * (-0.9) = 0.81
    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.81, state.value, 1e-15);
}

test "ForwardEuler converges to exact solution for exponential decay" {
    // dy/dt = -y, y(0) = 1 → y(t) = e^{-t}.
    // Forward Euler with small dt should approximate this.
    const Stepper = ForwardEuler(MockExplicitSystem);
    var state = MockExplicitSystem.State{ .value = 1.0 };

    const dt = 0.001;
    const steps = 1000;
    const t_final = @as(f64, @floatFromInt(steps)) * dt;

    for (0..steps) |_| {
        try Stepper.step(testing.allocator, &state, dt);
    }

    const exact = @exp(-t_final);
    // First-order method: error ~ O(dt) ~ 1e-3.
    try testing.expectApproxEqAbs(exact, state.value, 1e-2);
}

test "ForwardEuler works through TimeStepper wrapper" {
    const Stepper = time_stepper.TimeStepper(ForwardEuler(MockExplicitSystem));
    var state = MockExplicitSystem.State{ .value = 1.0 };

    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.9, state.value, 1e-15);
}

test "ForwardEuler on zero state is a no-op" {
    const Stepper = ForwardEuler(MockExplicitSystem);
    var state = MockExplicitSystem.State{ .value = 0.0 };

    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.0, state.value, 1e-15);
}

test "ForwardEuler and Leapfrog are distinct integrator types" {
    // This test validates that the interface isn't leapfrog-shaped:
    // ForwardEuler requires ExplicitSystem (forward_step), not
    // LeapfrogSystem (first_half_step + second_half_step).
    const leapfrog = @import("leapfrog.zig");

    // A type that satisfies ExplicitSystem but NOT LeapfrogSystem.
    const ExplicitOnly = struct {
        pub const State = struct { value: f64 };
        pub fn forward_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
            state.value += dt * (-state.value);
        }
    };

    // ForwardEuler accepts it.
    comptime ExplicitSystem(ExplicitOnly);
    const EulerStepper = ForwardEuler(ExplicitOnly);
    comptime time_stepper.TimeStepStrategy(EulerStepper);

    // Leapfrog rejects it (no first_half_step / second_half_step).
    // comptime leapfrog.LeapfrogSystem(ExplicitOnly);
    // expected: @compileError(...)
    _ = leapfrog;
}
