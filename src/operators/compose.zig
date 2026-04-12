//! Operator composition for type-checked chaining of DEC operators.
//!
//! The `chain` function applies a sequence of operators (d, ★, ★⁻¹, Δ, ...)
//! to a cochain, managing intermediate allocations automatically. Degree and
//! duality mismatches in the operator chain are compile errors.
//!
//! Example — the codifferential δ = ★⁻¹ ∘ d ∘ ★ on a primal 1-form:
//!
//!     const hodge_star = @import("hodge_star.zig");
//!     const exterior_derivative = @import("exterior_derivative.zig");
//!
//!     var result = try chain(allocator, omega, .{
//!         hodge_star.hodge_star,
//!         exterior_derivative.exterior_derivative,
//!         hodge_star.hodge_star_inverse,
//!     });
//!     defer result.deinit(allocator);

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const exterior_derivative = @import("exterior_derivative.zig");
const hodge_star = @import("hodge_star.zig");

/// Compute the result type of chaining `ops[start..]` starting from `InputType`.
///
/// Each operator in the tuple is a generic function `fn(Allocator, anytype) !U`;
/// we resolve each step's concrete return type by threading the output type of
/// each operator into the next, using `@TypeOf` on the call expression (which
/// resolves types without evaluating the function body).
fn ChainResultFrom(comptime InputType: type, comptime ops: anytype, comptime start: comptime_int) type {
    if (start >= ops.len) return InputType;
    const next = @typeInfo(@TypeOf(ops[start](
        @as(std.mem.Allocator, undefined),
        @as(InputType, undefined),
    ))).error_union.payload;
    return ChainResultFrom(next, ops, start + 1);
}

fn ChainResult(comptime InputType: type, comptime ops: anytype) type {
    return ChainResultFrom(InputType, ops, 0);
}

/// Apply a sequence of operators to a cochain, left to right.
///
/// Each operator must accept `(Allocator, previous_output)` and return
/// an error union of a cochain. Intermediate allocations are freed
/// automatically; only the final result is returned to the caller.
///
/// Degree or duality mismatches between adjacent operators in the chain
/// are compile errors — the types thread through at comptime.
pub fn chain(
    allocator: std.mem.Allocator,
    input: anytype,
    comptime ops: anytype,
) !ChainResult(@TypeOf(input), ops) {
    comptime {
        if (ops.len == 0) @compileError("chain requires at least one operator");
    }
    return chainFrom(allocator, input, ops, 0);
}

/// Index-based recursive implementation — applies ops[start], then chains the rest.
fn chainFrom(
    allocator: std.mem.Allocator,
    input: anytype,
    comptime ops: anytype,
    comptime start: comptime_int,
) !ChainResultFrom(@TypeOf(input), ops, start) {
    if (start >= ops.len) return input;

    if (start == ops.len - 1) {
        return ops[start](allocator, input);
    }

    // Apply the current operator.
    var intermediate = try ops[start](allocator, input);
    errdefer intermediate.deinit(allocator);

    // Recursively apply the remaining operators.
    const result = try chainFrom(allocator, intermediate, ops, start + 1);

    // Free the intermediate — the next operator has already read from it.
    intermediate.deinit(allocator);

    return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const PrimalC1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const PrimalC2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);

test "compile-time: chain type deduction for ★ ∘ d" {
    comptime {
        const DualC1 = cochain.Cochain(Mesh2D, 1, cochain.Dual);
        const Result = ChainResult(PrimalC0, .{ exterior_derivative.exterior_derivative, hodge_star.hodge_star });
        try testing.expect(Result == DualC1);
    }
}

test "compile-time: chain type deduction for ★⁻¹ ∘ d ∘ ★ (codifferential)" {
    comptime {
        // δ₁: primal 1-form → primal 0-form
        const Result = ChainResult(PrimalC1, .{ hodge_star.hodge_star, exterior_derivative.exterior_derivative, hodge_star.hodge_star_inverse });
        try testing.expect(Result == PrimalC0);
    }
}

test "chain with single operator equals direct call" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    const coords = mesh.vertices.slice().items(.coords);
    for (omega.values, coords) |*v, c| v.* = c[0];

    var direct = try exterior_derivative.exterior_derivative(allocator, omega);
    defer direct.deinit(allocator);

    var chained = try chain(allocator, omega, .{exterior_derivative.exterior_derivative});
    defer chained.deinit(allocator);

    for (direct.values, chained.values) |d, c| {
        try testing.expectApproxEqAbs(d, c, 1e-15);
    }
}

test "codifferential δ₁ via chain matches manual ★⁻¹ ∘ d ∘ ★" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_C047);

    for (0..100) |_| {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        // Manual: ★ → d → ★⁻¹
        var star_omega = try hodge_star.hodge_star(allocator, omega);
        defer star_omega.deinit(allocator);
        var d_star = try exterior_derivative.exterior_derivative(allocator, star_omega);
        defer d_star.deinit(allocator);
        var manual = try hodge_star.hodge_star_inverse(allocator, d_star);
        defer manual.deinit(allocator);

        // Chain
        var composed = try chain(allocator, omega, .{
            hodge_star.hodge_star,
            exterior_derivative.exterior_derivative,
            hodge_star.hodge_star_inverse,
        });
        defer composed.deinit(allocator);

        for (manual.values, composed.values) |m, c| {
            try testing.expectApproxEqAbs(m, c, 1e-15);
        }
    }
}

test "δd via chain matches laplacian on 0-forms" {
    // Δ₀ = δd = ★⁻¹ ∘ d ∘ ★ ∘ d (dδ vanishes for k=0).
    // Verify the chain produces the same result as the assembled context path.
    const context_mod = @import("context.zig");
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const operator_context = try context_mod.OperatorContext(Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.laplacian(0);

    var rng = std.Random.DefaultPrng.init(0xDEC_C049);

    for (0..100) |_| {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        // Via assembled operator context
        var lap = try (try operator_context.laplacian(0)).apply(allocator, omega);
        defer lap.deinit(allocator);

        // Via chain: d → ★ → d̃ → ★⁻¹ (which is ★⁻¹ ∘ d̃ ∘ ★ ∘ d = δd)
        var via_chain = try chain(allocator, omega, .{
            exterior_derivative.exterior_derivative,
            hodge_star.hodge_star,
            exterior_derivative.exterior_derivative,
            hodge_star.hodge_star_inverse,
        });
        defer via_chain.deinit(allocator);

        for (lap.values, via_chain.values) |l, c| {
            try testing.expectApproxEqAbs(l, c, 1e-12);
        }
    }
}
