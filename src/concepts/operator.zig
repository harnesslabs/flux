//! Comptime operator concept and generic wrapper.
//!
//! Two layers (following the TimeStepper pattern):
//!   1. `OperatorConcept` — validates that a type has the required `apply`
//!      declaration, plus comptime metadata (`InputType`, `OutputType`).
//!   2. Future: generic composition infrastructure that accepts any conforming
//!      operator type.
//!
//! An operator is a discrete linear map between cochain spaces:
//!   T: Ωᵏ(duality) → Ωˡ(duality')
//!
//! Conforming types must declare:
//!   - `pub const InputType: type`  — the cochain type consumed
//!   - `pub const OutputType: type` — the cochain type produced
//!   - `pub fn apply(std.mem.Allocator, InputType) !OutputType`

const std = @import("std");
const testing = std.testing;

// ═══════════════════════════════════════════════════════════════════════════
// OperatorConcept — comptime validation
// ═══════════════════════════════════════════════════════════════════════════

/// Validate that `Op` satisfies the Operator concept at compile time.
///
/// A conforming operator must declare:
///   - `pub const InputType: type`  — the cochain type this operator consumes
///   - `pub const OutputType: type` — the cochain type this operator produces
///   - `pub fn apply(std.mem.Allocator, InputType) !OutputType`
///
/// Produces a descriptive `@compileError` on violation.
pub fn OperatorConcept(comptime Op: type) void {
    _ = Op;
    @compileError("OperatorConcept not yet implemented");
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — OperatorConcept
// ═══════════════════════════════════════════════════════════════════════════

/// A minimal conforming operator for testing — identity on f64 slices.
const MockOperator = struct {
    pub const InputType = MockCochain;
    pub const OutputType = MockCochain;

    pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
        const values = try allocator.alloc(f64, input.values.len);
        @memcpy(values, input.values);
        return OutputType{ .values = values };
    }
};

const MockCochain = struct {
    values: []f64,

    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.values);
    }
};

test "OperatorConcept accepts a conforming type" {
    comptime OperatorConcept(MockOperator);
}

test "OperatorConcept accepts an operator with extra declarations" {
    const ExtendedOp = struct {
        pub const InputType = MockCochain;
        pub const OutputType = MockCochain;

        pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
            const values = try allocator.alloc(f64, input.values.len);
            @memcpy(values, input.values);
            return OutputType{ .values = values };
        }

        /// Extra method — concept should not reject this.
        pub fn name() []const u8 {
            return "extended";
        }
    };
    comptime OperatorConcept(ExtendedOp);
}

// ── Negative tests (compile-time rejection) ──────────────────────────────
//
// These verify that non-conforming types produce compile errors. Since
// @compileError is fatal, we cannot test them in CI. Verification:
// uncomment the `comptime` call and confirm it fails.

test "OperatorConcept rejects type missing InputType" {
    const NoInput = struct {
        pub const OutputType = MockCochain;
        pub fn apply(_: std.mem.Allocator, _: MockCochain) !MockCochain {
            unreachable;
        }
    };

    // comptime OperatorConcept(NoInput);
    // expected: @compileError("OperatorConcept requires a 'pub const InputType: type' declaration ...")
    _ = NoInput;
}

test "OperatorConcept rejects type missing OutputType" {
    const NoOutput = struct {
        pub const InputType = MockCochain;
        pub fn apply(_: std.mem.Allocator, _: MockCochain) !MockCochain {
            unreachable;
        }
    };

    // comptime OperatorConcept(NoOutput);
    // expected: @compileError("OperatorConcept requires a 'pub const OutputType: type' declaration ...")
    _ = NoOutput;
}

test "OperatorConcept rejects type missing apply function" {
    const NoApply = struct {
        pub const InputType = MockCochain;
        pub const OutputType = MockCochain;
    };

    // comptime OperatorConcept(NoApply);
    // expected: @compileError("OperatorConcept requires a 'pub fn apply(...)' declaration ...")
    _ = NoApply;
}

test "OperatorConcept rejects apply with wrong parameter count" {
    const WrongArity = struct {
        pub const InputType = MockCochain;
        pub const OutputType = MockCochain;

        pub fn apply(_: MockCochain) !MockCochain {
            unreachable;
        }
    };

    // comptime OperatorConcept(WrongArity);
    // expected: @compileError("OperatorConcept: 'apply' must take exactly 2 parameters ...")
    _ = WrongArity;
}

test "OperatorConcept rejects apply with wrong input type" {
    const WrongInput = struct {
        pub const InputType = MockCochain;
        pub const OutputType = MockCochain;

        pub fn apply(_: std.mem.Allocator, _: u32) !MockCochain {
            unreachable;
        }
    };

    // comptime OperatorConcept(WrongInput);
    // expected: @compileError("OperatorConcept: 'apply' parameter 1 must be InputType ...")
    _ = WrongInput;
}

test "OperatorConcept rejects apply with wrong return type" {
    const WrongReturn = struct {
        pub const InputType = MockCochain;
        pub const OutputType = MockCochain;

        pub fn apply(_: std.mem.Allocator, _: MockCochain) !u32 {
            unreachable;
        }
    };

    // comptime OperatorConcept(WrongReturn);
    // expected: @compileError("OperatorConcept: 'apply' must return !OutputType ...")
    _ = WrongReturn;
}

test "MockOperator.apply copies input values" {
    const allocator = testing.allocator;

    var input_buf = [_]f64{ 1.0, 2.0, 3.0 };
    const input = MockCochain{ .values = &input_buf };

    var output = try MockOperator.apply(allocator, input);
    defer allocator.free(output.values);

    try testing.expectEqual(@as(usize, 3), output.values.len);
    try testing.expectApproxEqAbs(1.0, output.values[0], 1e-15);
    try testing.expectApproxEqAbs(2.0, output.values[1], 1e-15);
    try testing.expectApproxEqAbs(3.0, output.values[2], 1e-15);
}
