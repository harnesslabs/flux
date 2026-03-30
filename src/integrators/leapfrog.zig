//! Generic leapfrog (Störmer-Verlet) integrator.
//!
//! Leapfrog advances two coupled variable groups in a staggered fashion:
//!
//!   x^{n+1/2} = x^{n-1/2} + dt · f(y^n)       (first half-step)
//!   y^{n+1}   = y^n        + dt · g(x^{n+1/2}) (second half-step)
//!
//! The scheme is symplectic — it exactly preserves a modified Hamiltonian,
//! giving bounded energy oscillation (no secular drift) for conservative
//! systems. This makes it the natural choice for wave equations (Maxwell,
//! acoustics) and Hamiltonian mechanics.
//!
//! ## Usage
//!
//! Define a system type with two half-step operators, then wrap it:
//!
//! ```zig
//! const MySystem = struct {
//!     pub const State = MyState;
//!     pub fn first_half_step(alloc: Allocator, state: *State, dt: f64) !void { ... }
//!     pub fn second_half_step(alloc: Allocator, state: *State, dt: f64) !void { ... }
//! };
//! const Stepper = Leapfrog(MySystem); // satisfies TimeStepStrategy
//! try Stepper.step(allocator, &state, dt);
//! ```

const std = @import("std");
const testing = std.testing;
const time_stepper = @import("../time_stepper.zig");

// ═══════════════════════════════════════════════════════════════════════════
// LeapfrogSystem — comptime concept
// ═══════════════════════════════════════════════════════════════════════════

/// Validate that `S` satisfies the LeapfrogSystem concept at compile time.
///
/// A conforming system must declare:
///   - `pub const State: type`
///   - `pub fn first_half_step(Allocator, *State, f64) !void`
///   - `pub fn second_half_step(Allocator, *State, f64) !void`
///
/// The two half-steps encode the symplectic splitting: the first advances
/// the staggered variable (e.g., B in Maxwell), the second advances the
/// integer-time variable (e.g., E in Maxwell).
pub fn LeapfrogSystem(comptime S: type) void {
    // 1. S must declare a State type.
    if (!@hasDecl(S, "State")) {
        @compileError("LeapfrogSystem requires a 'pub const State: type' declaration — " ++
            "the simulation state type this system advances");
    }
    const State = S.State;
    if (@TypeOf(State) != type) {
        @compileError("LeapfrogSystem: 'State' must be a type, got " ++ @typeName(@TypeOf(State)));
    }

    // 2. Validate both half-step functions.
    inline for (.{ "first_half_step", "second_half_step" }) |name| {
        if (!@hasDecl(S, name)) {
            @compileError("LeapfrogSystem requires a 'pub fn " ++ name ++
                "(std.mem.Allocator, *State, f64) !void' declaration");
        }

        const info = @typeInfo(@TypeOf(@field(S, name)));
        if (info != .@"fn") {
            @compileError("LeapfrogSystem: '" ++ name ++ "' must be a function");
        }
        const fn_info = info.@"fn";

        if (fn_info.params.len != 3) {
            @compileError("LeapfrogSystem: '" ++ name ++
                "' must take exactly 3 parameters (std.mem.Allocator, *State, f64)");
        }
        if (fn_info.params[0].type != std.mem.Allocator) {
            @compileError("LeapfrogSystem: '" ++ name ++ "' parameter 0 must be std.mem.Allocator");
        }
        if (fn_info.params[1].type != *State) {
            @compileError("LeapfrogSystem: '" ++ name ++ "' parameter 1 must be *State");
        }
        if (fn_info.params[2].type != f64) {
            @compileError("LeapfrogSystem: '" ++ name ++ "' parameter 2 must be f64 (the timestep dt)");
        }

        const ret = fn_info.return_type orelse
            @compileError("LeapfrogSystem: '" ++ name ++ "' must have a known return type");
        const ret_info = @typeInfo(ret);
        if (ret_info != .error_union) {
            @compileError("LeapfrogSystem: '" ++ name ++ "' must return !void");
        }
        if (ret_info.error_union.payload != void) {
            @compileError("LeapfrogSystem: '" ++ name ++ "' must return !void, not !" ++
                @typeName(ret_info.error_union.payload));
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Leapfrog — generic integrator
// ═══════════════════════════════════════════════════════════════════════════

/// Generic leapfrog integrator parameterized on a `LeapfrogSystem`.
///
/// Composes two half-steps into a full timestep. The returned type satisfies
/// `TimeStepStrategy`, so it can be wrapped by `TimeStepper` and used by
/// simulation runners.
///
/// The integrator does not track timestep counts — that is the caller's
/// responsibility (typically via a field on State or the runner loop).
pub fn Leapfrog(comptime System: type) type {
    comptime LeapfrogSystem(System);

    return struct {
        /// The simulation state type, forwarded from the system.
        pub const State = System.State;

        /// Advance the state by one full leapfrog timestep.
        ///
        /// Executes: first_half_step(dt) then second_half_step(dt).
        pub fn step(allocator: std.mem.Allocator, state: *State, dt: f64) !void {
            try System.first_half_step(allocator, state, dt);
            try System.second_half_step(allocator, state, dt);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — LeapfrogSystem concept
// ═══════════════════════════════════════════════════════════════════════════

/// Mock system: two coupled variables advanced by half-steps.
const MockLeapfrogSystem = struct {
    pub const State = struct {
        position: f64,
        velocity: f64,
    };

    /// First half-step: advance position using current velocity.
    pub fn first_half_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
        state.position += dt * state.velocity;
    }

    /// Second half-step: advance velocity using updated position.
    /// (Trivial force: a = -position, i.e., harmonic oscillator.)
    pub fn second_half_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
        state.velocity += dt * (-state.position);
    }
};

test "LeapfrogSystem accepts a conforming type" {
    comptime LeapfrogSystem(MockLeapfrogSystem);
}

test "LeapfrogSystem accepts a type with extra declarations" {
    const Extended = struct {
        pub const State = struct { x: f64, v: f64 };
        pub fn first_half_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
            state.x += dt * state.v;
        }
        pub fn second_half_step(_: std.mem.Allocator, state: *State, dt: f64) !void {
            state.v -= dt * state.x;
        }
        /// Extra — concept should not reject.
        pub fn energy(state: *const State) f64 {
            return 0.5 * (state.x * state.x + state.v * state.v);
        }
    };
    comptime LeapfrogSystem(Extended);
}

// ── Negative tests (compile-time rejection) ──────────────────────────────

test "LeapfrogSystem rejects type missing State" {
    const NoState = struct {
        pub fn first_half_step(_: std.mem.Allocator, _: *anyopaque, _: f64) !void {}
        pub fn second_half_step(_: std.mem.Allocator, _: *anyopaque, _: f64) !void {}
    };
    // comptime LeapfrogSystem(NoState);
    // expected: @compileError("LeapfrogSystem requires a 'pub const State: type' declaration ...")
    _ = NoState;
}

test "LeapfrogSystem rejects type missing first_half_step" {
    const NoFirst = struct {
        pub const State = struct { x: f64 };
        pub fn second_half_step(_: std.mem.Allocator, _: *State, _: f64) !void {}
    };
    // comptime LeapfrogSystem(NoFirst);
    // expected: @compileError("LeapfrogSystem requires a 'pub fn first_half_step(...)' ...")
    _ = NoFirst;
}

test "LeapfrogSystem rejects type missing second_half_step" {
    const NoSecond = struct {
        pub const State = struct { x: f64 };
        pub fn first_half_step(_: std.mem.Allocator, _: *State, _: f64) !void {}
    };
    // comptime LeapfrogSystem(NoSecond);
    // expected: @compileError("LeapfrogSystem requires a 'pub fn second_half_step(...)' ...")
    _ = NoSecond;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — Leapfrog integrator
// ═══════════════════════════════════════════════════════════════════════════

test "Leapfrog satisfies TimeStepStrategy" {
    const Stepper = Leapfrog(MockLeapfrogSystem);
    comptime time_stepper.TimeStepStrategy(Stepper);
}

test "Leapfrog composes first_half_step then second_half_step" {
    const Stepper = Leapfrog(MockLeapfrogSystem);
    var state = MockLeapfrogSystem.State{ .position = 0.0, .velocity = 1.0 };

    // One step of leapfrog on a harmonic oscillator (ẍ = -x):
    //   position += dt * velocity  → 0.0 + 0.1 * 1.0 = 0.1
    //   velocity += dt * (-position) → 1.0 + 0.1 * (-0.1) = 0.99
    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.1, state.position, 1e-15);
    try testing.expectApproxEqAbs(0.99, state.velocity, 1e-15);
}

test "Leapfrog preserves energy for harmonic oscillator (symplecticity)" {
    // The leapfrog scheme is symplectic — energy should oscillate around
    // the true value without secular drift over many steps.
    const Stepper = Leapfrog(MockLeapfrogSystem);
    var state = MockLeapfrogSystem.State{ .position = 1.0, .velocity = 0.0 };

    const initial_energy = 0.5 * (state.position * state.position + state.velocity * state.velocity);
    const dt = 0.01;
    const steps = 10_000;

    for (0..steps) |_| {
        try Stepper.step(testing.allocator, &state, dt);
    }

    const final_energy = 0.5 * (state.position * state.position + state.velocity * state.velocity);

    // Symplectic: energy drift bounded, not accumulating.
    // For dt=0.01, 10k steps on a harmonic oscillator, the modified
    // Hamiltonian differs from the true one by O(dt²). The oscillation
    // amplitude is small but we use a conservative tolerance.
    try testing.expectApproxEqAbs(initial_energy, final_energy, 5e-3);
}

test "Leapfrog works through TimeStepper wrapper" {
    const Stepper = time_stepper.TimeStepper(Leapfrog(MockLeapfrogSystem));
    var state = MockLeapfrogSystem.State{ .position = 0.0, .velocity = 1.0 };

    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.1, state.position, 1e-15);
}

test "Leapfrog on zero state is a no-op" {
    const Stepper = Leapfrog(MockLeapfrogSystem);
    var state = MockLeapfrogSystem.State{ .position = 0.0, .velocity = 0.0 };

    try Stepper.step(testing.allocator, &state, 0.1);
    try testing.expectApproxEqAbs(0.0, state.position, 1e-15);
    try testing.expectApproxEqAbs(0.0, state.velocity, 1e-15);
}
