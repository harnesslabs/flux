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
    // 1. Op must declare an InputType type.
    if (!@hasDecl(Op, "InputType")) {
        @compileError("OperatorConcept requires a 'pub const InputType: type' declaration — " ++
            "the cochain type this operator consumes");
    }
    const Input = Op.InputType;
    if (@TypeOf(Input) != type) {
        @compileError("OperatorConcept: 'InputType' must be a type, got " ++ @typeName(@TypeOf(Input)));
    }

    // 2. Op must declare an OutputType type.
    if (!@hasDecl(Op, "OutputType")) {
        @compileError("OperatorConcept requires a 'pub const OutputType: type' declaration — " ++
            "the cochain type this operator produces");
    }
    const Output = Op.OutputType;
    if (@TypeOf(Output) != type) {
        @compileError("OperatorConcept: 'OutputType' must be a type, got " ++ @typeName(@TypeOf(Output)));
    }

    // 3. Op must declare a pub apply function: fn(Allocator, InputType) !OutputType
    if (!@hasDecl(Op, "apply")) {
        @compileError("OperatorConcept requires a 'pub fn apply(std.mem.Allocator, InputType) !OutputType' " ++
            "declaration — apply the operator to an input cochain");
    }

    const apply_info = @typeInfo(@TypeOf(Op.apply));
    if (apply_info != .@"fn") {
        @compileError("OperatorConcept: 'apply' must be a function");
    }
    const fn_info = apply_info.@"fn";

    // Check parameter count and types.
    if (fn_info.params.len != 2) {
        @compileError("OperatorConcept: 'apply' must take exactly 2 parameters " ++
            "(std.mem.Allocator, InputType)");
    }
    if (fn_info.params[0].type != std.mem.Allocator) {
        @compileError("OperatorConcept: 'apply' parameter 0 must be std.mem.Allocator");
    }
    if (fn_info.params[1].type != Input) {
        @compileError("OperatorConcept: 'apply' parameter 1 must be InputType (" ++
            @typeName(Input) ++ ")");
    }

    // Check return type is an error union with payload OutputType.
    const ret = fn_info.return_type orelse
        @compileError("OperatorConcept: 'apply' must have a known return type");
    const ret_info = @typeInfo(ret);
    if (ret_info != .error_union) {
        @compileError("OperatorConcept: 'apply' must return !OutputType (an error union)");
    }
    if (ret_info.error_union.payload != Output) {
        @compileError("OperatorConcept: 'apply' must return !OutputType (!" ++
            @typeName(Output) ++ "), got !" ++ @typeName(ret_info.error_union.payload));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Operator wrappers — retrofit free functions into concept-conforming types
// ═══════════════════════════════════════════════════════════════════════════

/// Exterior derivative dₖ as a concept-conforming operator type.
///
/// Wraps the free function `exterior_derivative(allocator, input)` into a
/// type with `InputType`, `OutputType`, and `apply`. Parameterized on mesh
/// type, degree k, and duality so the cochain types are fully resolved at
/// comptime.
///
/// ```zig
/// const D0 = ExteriorDerivativeOperator(Mesh(2), 0, Primal);
/// comptime OperatorConcept(D0);  // passes
/// var result = try D0.apply(allocator, omega_0);
/// ```
pub fn ExteriorDerivativeOperator(comptime MeshType: type, comptime k: comptime_int, comptime D: type) type {
    const forms = @import("../forms/cochain.zig");
    const ext = @import("../operators/exterior_derivative.zig");
    const Input = forms.Cochain(MeshType, k, D);
    const Output = forms.Cochain(MeshType, k + 1, D);

    const Op = struct {
        pub const InputType = Input;
        pub const OutputType = Output;

        pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
            return ext.exterior_derivative(allocator, input);
        }
    };

    // Gate: the constructed type must satisfy the concept.
    comptime OperatorConcept(Op);

    return Op;
}

/// Hodge star ★ₖ as a concept-conforming operator type.
///
/// Maps primal k-cochains to dual (n−k)-cochains. Wraps `hodge_star()`.
pub fn HodgeStarOperator(comptime MeshType: type, comptime k: comptime_int) type {
    const forms = @import("../forms/cochain.zig");
    const hs = @import("../operators/hodge_star.zig");
    const n = MeshType.topological_dimension;
    const Input = forms.Cochain(MeshType, k, forms.Primal);
    const Output = forms.Cochain(MeshType, n - k, forms.Dual);

    const Op = struct {
        pub const InputType = Input;
        pub const OutputType = Output;

        pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
            return hs.hodge_star(allocator, input);
        }
    };

    comptime OperatorConcept(Op);

    return Op;
}

/// Inverse Hodge star ★⁻¹ as a concept-conforming operator type.
///
/// Maps dual (n−k)-cochains to primal k-cochains. Wraps `hodge_star_inverse()`.
pub fn HodgeStarInverseOperator(comptime MeshType: type, comptime k: comptime_int) type {
    const forms = @import("../forms/cochain.zig");
    const hs = @import("../operators/hodge_star.zig");
    const n = MeshType.topological_dimension;
    const Input = forms.Cochain(MeshType, n - k, forms.Dual);
    const Output = forms.Cochain(MeshType, k, forms.Primal);

    const Op = struct {
        pub const InputType = Input;
        pub const OutputType = Output;

        pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
            return hs.hodge_star_inverse(allocator, input);
        }
    };

    comptime OperatorConcept(Op);

    return Op;
}

/// Hodge Laplacian Δₖ as a concept-conforming operator type.
///
/// Maps primal k-cochains to primal k-cochains. Wraps `laplacian()`.
pub fn LaplacianOperator(comptime MeshType: type, comptime k: comptime_int) type {
    const forms = @import("../forms/cochain.zig");
    const lap = @import("../operators/laplacian.zig");
    const CochainType = forms.Cochain(MeshType, k, forms.Primal);

    const Op = struct {
        pub const InputType = CochainType;
        pub const OutputType = CochainType;

        pub fn apply(allocator: std.mem.Allocator, input: InputType) !OutputType {
            return lap.laplacian(allocator, input);
        }
    };

    comptime OperatorConcept(Op);

    return Op;
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

// ── Operator wrapper tests — concept conformance ────────────────────────

const topology = @import("../topology/mesh.zig");
const Mesh2D = topology.Mesh(2);

test "ExteriorDerivativeOperator(Mesh(2), 0, Primal) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const D0 = ExteriorDerivativeOperator(Mesh2D, 0, forms.Primal);
    comptime OperatorConcept(D0);

    // Verify the type metadata is correct.
    comptime {
        std.debug.assert(D0.InputType.degree == 0);
        std.debug.assert(D0.OutputType.degree == 1);
        std.debug.assert(D0.InputType.duality == forms.Primal);
        std.debug.assert(D0.OutputType.duality == forms.Primal);
    }
}

test "ExteriorDerivativeOperator(Mesh(2), 1, Primal) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const D1 = ExteriorDerivativeOperator(Mesh2D, 1, forms.Primal);
    comptime OperatorConcept(D1);

    comptime {
        std.debug.assert(D1.InputType.degree == 1);
        std.debug.assert(D1.OutputType.degree == 2);
    }
}

test "ExteriorDerivativeOperator(Mesh(2), 0, Dual) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const D0Dual = ExteriorDerivativeOperator(Mesh2D, 0, forms.Dual);
    comptime OperatorConcept(D0Dual);

    comptime {
        std.debug.assert(D0Dual.InputType.duality == forms.Dual);
        std.debug.assert(D0Dual.OutputType.duality == forms.Dual);
    }
}

test "HodgeStarOperator(Mesh(2), 0) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const Star0 = HodgeStarOperator(Mesh2D, 0);
    comptime OperatorConcept(Star0);

    comptime {
        std.debug.assert(Star0.InputType.degree == 0);
        std.debug.assert(Star0.InputType.duality == forms.Primal);
        std.debug.assert(Star0.OutputType.degree == 2);
        std.debug.assert(Star0.OutputType.duality == forms.Dual);
    }
}

test "HodgeStarOperator(Mesh(2), 1) satisfies OperatorConcept" {
    const Star1 = HodgeStarOperator(Mesh2D, 1);
    comptime OperatorConcept(Star1);
}

test "HodgeStarOperator(Mesh(2), 2) satisfies OperatorConcept" {
    const Star2 = HodgeStarOperator(Mesh2D, 2);
    comptime OperatorConcept(Star2);
}

test "HodgeStarInverseOperator(Mesh(2), 0) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const StarInv0 = HodgeStarInverseOperator(Mesh2D, 0);
    comptime OperatorConcept(StarInv0);

    comptime {
        // ★⁻¹₀: dual 2-form → primal 0-form
        std.debug.assert(StarInv0.InputType.degree == 2);
        std.debug.assert(StarInv0.InputType.duality == forms.Dual);
        std.debug.assert(StarInv0.OutputType.degree == 0);
        std.debug.assert(StarInv0.OutputType.duality == forms.Primal);
    }
}

test "HodgeStarInverseOperator(Mesh(2), 1) satisfies OperatorConcept" {
    const StarInv1 = HodgeStarInverseOperator(Mesh2D, 1);
    comptime OperatorConcept(StarInv1);
}

test "LaplacianOperator(Mesh(2), 0) satisfies OperatorConcept" {
    const forms = @import("../forms/cochain.zig");
    const Lap0 = LaplacianOperator(Mesh2D, 0);
    comptime OperatorConcept(Lap0);

    comptime {
        // Δ₀: primal 0-form → primal 0-form (same type in and out)
        std.debug.assert(Lap0.InputType == Lap0.OutputType);
        std.debug.assert(Lap0.InputType.degree == 0);
        std.debug.assert(Lap0.InputType.duality == forms.Primal);
    }
}

test "LaplacianOperator(Mesh(2), 1) satisfies OperatorConcept" {
    const Lap1 = LaplacianOperator(Mesh2D, 1);
    comptime OperatorConcept(Lap1);
}

test "LaplacianOperator(Mesh(2), 2) satisfies OperatorConcept" {
    const Lap2 = LaplacianOperator(Mesh2D, 2);
    comptime OperatorConcept(Lap2);
}

// ── Runtime test — verify wrapper delegates to real implementation ───────

test "ExteriorDerivativeOperator.apply produces correct result on real mesh" {
    const forms = @import("../forms/cochain.zig");
    const allocator = testing.allocator;

    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // Constant 0-form: d(constant) = 0
    const Cochain0 = forms.Cochain(Mesh2D, 0, forms.Primal);
    var omega = try Cochain0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    @memset(omega.values, 1.0);

    const D0 = ExteriorDerivativeOperator(Mesh2D, 0, forms.Primal);
    var result = try D0.apply(allocator, omega);
    defer result.deinit(allocator);

    // d(constant) = 0 exactly
    for (result.values) |v| {
        try testing.expectApproxEqAbs(0.0, v, 1e-15);
    }
}

test "MockOperator.apply copies input values" {
    const allocator = testing.allocator;

    var input_buf = [_]f64{ 1.0, 2.0, 3.0 };
    const input = MockCochain{ .values = &input_buf };

    const output = try MockOperator.apply(allocator, input);
    defer allocator.free(output.values);

    try testing.expectEqual(@as(usize, 3), output.values.len);
    try testing.expectApproxEqAbs(1.0, output.values[0], 1e-15);
    try testing.expectApproxEqAbs(2.0, output.values[1], 1e-15);
    try testing.expectApproxEqAbs(3.0, output.values[2], 1e-15);
}
