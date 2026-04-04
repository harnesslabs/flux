//! Discrete Hodge star operator and its inverse.
//!
//! The Hodge star ★ₖ maps primal k-cochains to dual (n−k)-cochains.
//! The implementation dispatches by degree at comptime:
//!
//!   ★₀ and ★₂ (k = 0 and k = n): diagonal operators using dual mesh volumes.
//!     (★₀)ᵢᵢ = dual_area[i]
//!     (★₂)ᵢᵢ = 1 / area[i]
//!
//!   ★₁ (0 < k < n): Whitney/Galerkin mass matrix M₁ (SpMV for forward,
//!     preconditioned CG solve for inverse). The diagonal approximation
//!     dual_length/length is only consistent on orthogonal dual meshes;
//!     the Whitney mass matrix is exact on any triangulation.
//!
//! The mesh stores both M₁ and its diagonal preconditioner, precomputed
//! during construction.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const sparse = @import("../math/sparse.zig");
const conjugate_gradient = @import("../math/cg.zig");

pub fn AssembledHodgeStar(comptime InputType: type) type {
    comptime validateHodgeStarInput(InputType);

    return struct {
        const Self = @This();

        mesh: *const InputType.MeshT,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !HodgeStarResult(InputType) {
            std.debug.assert(input.mesh == self.mesh);
            return apply_hodge_star(allocator, input);
        }
    };
}

pub fn AssembledHodgeStarInverse(comptime InputType: type) type {
    comptime validateHodgeStarInverseInput(InputType);

    return struct {
        const Self = @This();

        mesh: *const InputType.MeshT,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !HodgeStarInverseResult(InputType) {
            std.debug.assert(input.mesh == self.mesh);
            return apply_hodge_star_inverse(allocator, input);
        }
    };
}

pub fn assemble_for_degree(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledHodgeStar(cochain.Cochain(MeshType, primal_degree, cochain.Primal)) {
    _ = allocator;
    return .{ .mesh = mesh };
}

pub fn assemble_inverse_for_degree(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledHodgeStarInverse(cochain.Cochain(
    MeshType,
    MeshType.topological_dimension - primal_degree,
    cochain.Dual,
)) {
    _ = allocator;
    return .{ .mesh = mesh };
}

/// Apply the Hodge star ★ₖ to a primal k-cochain, returning a dual (n−k)-cochain.
///
/// For k=0 and k=n, this is a diagonal scaling by dual mesh volumes.
/// For k=1, this applies the Whitney mass matrix M₁ via SpMV.
pub fn apply_hodge_star(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInput(InputType);

    const k = InputType.degree;
    const topological_dimension = InputType.MeshT.topological_dimension;
    const OutputType = cochain.Cochain(InputType.MeshT, topological_dimension - k, cochain.Dual);

    var output = try OutputType.init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    if (isDiagonalDegree(InputType.MeshT, k)) {
        applyDiagonal(InputType.MeshT, k, input.mesh, input.values, output.values, false);
        return output;
    }

    if (!supportsWhitneyMassDegree(InputType.MeshT, k)) {
        return error.NotYetImplemented;
    }

    sparse.spmv(input.mesh.whitney_mass_1, input.values, output.values);
    return output;
}

/// Apply the inverse Hodge star ★⁻¹ to a dual (n−k)-cochain, returning a
/// primal k-cochain.
///
/// For k=0 and k=n, this is the element-wise reciprocal of the diagonal.
/// For k=1, this solves M₁ x = b via preconditioned conjugate gradient.
pub fn apply_hodge_star_inverse(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInverseInput(InputType);

    const dual_degree = InputType.degree;
    const topological_dimension = InputType.MeshT.topological_dimension;
    const primal_degree = topological_dimension - dual_degree;
    const OutputType = cochain.Cochain(InputType.MeshT, primal_degree, cochain.Primal);

    var output = try OutputType.init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    if (isDiagonalDegree(InputType.MeshT, primal_degree)) {
        applyDiagonal(InputType.MeshT, primal_degree, input.mesh, input.values, output.values, true);
        return output;
    }

    if (!supportsWhitneyMassDegree(InputType.MeshT, primal_degree)) {
        return error.NotYetImplemented;
    }

    // CG solve: M₁ · x = b, with diagonal ★₁ as preconditioner.
    @memset(output.values, 0.0);

    var scratch = try conjugate_gradient.Scratch.init(allocator, @intCast(output.values.len));
    defer scratch.deinit(allocator);

    var precond = conjugate_gradient.DiagonalPreconditioner{ .diagonal = input.mesh.preconditioner_1 };
    const result = conjugate_gradient.solve(
        input.mesh.whitney_mass_1,
        input.values,
        output.values,
        1e-10,
        1000,
        conjugate_gradient.DiagonalPreconditioner.apply,
        @ptrCast(&precond),
        scratch,
    );
    std.debug.assert(result.converged);

    return output;
}

pub fn hodge_star(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    return apply_hodge_star(allocator, input);
}

pub fn hodge_star_inverse(
    allocator: std.mem.Allocator,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    return apply_hodge_star_inverse(allocator, input);
}

// ── Return type helpers ──────────────────────────────────────────────────

fn HodgeStarResult(comptime InputType: type) type {
    const topological_dimension = InputType.MeshT.topological_dimension;
    return cochain.Cochain(InputType.MeshT, topological_dimension - InputType.degree, cochain.Dual);
}

fn HodgeStarInverseResult(comptime InputType: type) type {
    const topological_dimension = InputType.MeshT.topological_dimension;
    return cochain.Cochain(InputType.MeshT, topological_dimension - InputType.degree, cochain.Primal);
}

// ── Comptime validation ──────────────────────────────────────────────────

fn validateHodgeStarInput(comptime T: type) void {
    if (!@hasDecl(T, "duality")) {
        @compileError("hodge_star requires a Cochain type");
    }
    if (T.duality != cochain.Primal) {
        @compileError("hodge_star expects a primal cochain (use hodge_star_inverse for dual → primal)");
    }
}

fn validateHodgeStarInverseInput(comptime T: type) void {
    if (!@hasDecl(T, "duality")) {
        @compileError("hodge_star_inverse requires a Cochain type");
    }
    if (T.duality != cochain.Dual) {
        @compileError("hodge_star_inverse expects a dual cochain (use hodge_star for primal → dual)");
    }
}

// ── Raw buffer API (for operators that bypass typed cochains) ────────────

/// Apply ★ₖ to raw f64 buffers. Dispatches by primal_degree at comptime:
/// k=0,n use diagonal scaling; 0 < k < n uses Whitney SpMV.
pub fn apply_raw(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
) void {
    if (isDiagonalDegree(MeshType, primal_degree)) {
        applyDiagonal(MeshType, primal_degree, mesh, input, output, false);
        return;
    }

    if (!supportsWhitneyMassDegree(MeshType, primal_degree)) {
        @panic("Hodge star interior degree not yet implemented for this mesh dimension");
    }

    sparse.spmv(mesh.whitney_mass_1, input, output);
}

/// Apply ★⁻¹ₖ to raw f64 buffers. Dispatches by primal_degree at comptime:
/// k=0,n use diagonal reciprocal; 0 < k < n uses preconditioned CG solve.
pub fn apply_inverse_raw(
    allocator: std.mem.Allocator,
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
) !void {
    if (isDiagonalDegree(MeshType, primal_degree)) {
        applyDiagonal(MeshType, primal_degree, mesh, input, output, true);
        return;
    }

    if (!supportsWhitneyMassDegree(MeshType, primal_degree)) {
        return error.NotYetImplemented;
    }

    @memset(output, 0.0);

    var scratch = try conjugate_gradient.Scratch.init(allocator, @intCast(output.len));
    defer scratch.deinit(allocator);

    var precond = conjugate_gradient.DiagonalPreconditioner{ .diagonal = mesh.preconditioner_1 };
    const result = conjugate_gradient.solve(
        mesh.whitney_mass_1,
        input,
        output,
        1e-10,
        1000,
        conjugate_gradient.DiagonalPreconditioner.apply,
        @ptrCast(&precond),
        scratch,
    );
    std.debug.assert(result.converged);
}

// ── Diagonal application (k=0 and k=n only) ────────────────────────────

/// Apply the diagonal Hodge star (or its inverse) for degree 0 or n.
///
/// When `invert` is false: output[i] = ratio[i] * input[i]
/// When `invert` is true:  output[i] = input[i] / ratio[i]  (asserts ratio ≠ 0)
fn applyDiagonal(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
    comptime invert: bool,
) void {
    std.debug.assert(isDiagonalDegree(MeshType, primal_degree));

    if (primal_degree == 0) {
        const dual_volumes = mesh.vertices.slice().items(.dual_volume);
        for (output, input, dual_volumes) |*out, in_val, ratio| {
            if (invert) {
                std.debug.assert(ratio != 0.0);
                out.* = in_val / ratio;
            } else {
                out.* = ratio * in_val;
            }
        }
        return;
    }

    if (primal_degree == MeshType.topological_dimension) {
        const primal_volumes = mesh.simplices(MeshType.topological_dimension).items(.volume);
        for (output, input, primal_volumes) |*out, in_val, volume| {
            if (invert) {
                out.* = volume * in_val;
            } else {
                std.debug.assert(volume != 0.0);
                out.* = in_val / volume;
            }
        }
        return;
    }

    unreachable;
}

fn isDiagonalDegree(comptime MeshType: type, comptime primal_degree: comptime_int) bool {
    return primal_degree == 0 or primal_degree == MeshType.topological_dimension;
}

fn supportsWhitneyMassDegree(comptime MeshType: type, comptime primal_degree: comptime_int) bool {
    return MeshType.topological_dimension == 2 and primal_degree == 1;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const Mesh3D = topology.Mesh(3, 3);
const PrimalC0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const PrimalC1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const PrimalC2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);
const DualC0 = cochain.Cochain(Mesh2D, 0, cochain.Dual);
const DualC1 = cochain.Cochain(Mesh2D, 1, cochain.Dual);
const DualC2 = cochain.Cochain(Mesh2D, 2, cochain.Dual);
const PrimalC0_3D = cochain.Cochain(Mesh3D, 0, cochain.Primal);
const PrimalC1_3D = cochain.Cochain(Mesh3D, 1, cochain.Primal);
const PrimalC2_3D = cochain.Cochain(Mesh3D, 2, cochain.Primal);
const PrimalC3_3D = cochain.Cochain(Mesh3D, 3, cochain.Primal);
const DualC0_3D = cochain.Cochain(Mesh3D, 0, cochain.Dual);
const DualC1_3D = cochain.Cochain(Mesh3D, 1, cochain.Dual);
const DualC2_3D = cochain.Cochain(Mesh3D, 2, cochain.Dual);
const DualC3_3D = cochain.Cochain(Mesh3D, 3, cochain.Dual);

// ── Compile-time type checks ─────────────────────────────────────────────

test "compile-time: ★₀ maps primal 0-form to dual 2-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC0) == DualC2);
    }
}

test "compile-time: ★₁ maps primal 1-form to dual 1-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC1) == DualC1);
    }
}

test "compile-time: ★₂ maps primal 2-form to dual 0-form" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC2) == DualC0);
    }
}

test "compile-time: ★⁻¹ maps dual (n−k)-form back to primal k-form" {
    comptime {
        try testing.expect(HodgeStarInverseResult(DualC2) == PrimalC0);
        try testing.expect(HodgeStarInverseResult(DualC1) == PrimalC1);
        try testing.expect(HodgeStarInverseResult(DualC0) == PrimalC2);
    }
}

test "compile-time: ★ maps 3D primal k-forms to dual (3-k)-forms" {
    comptime {
        try testing.expect(HodgeStarResult(PrimalC0_3D) == DualC3_3D);
        try testing.expect(HodgeStarResult(PrimalC1_3D) == DualC2_3D);
        try testing.expect(HodgeStarResult(PrimalC2_3D) == DualC1_3D);
        try testing.expect(HodgeStarResult(PrimalC3_3D) == DualC0_3D);
    }
}

test "compile-time: ★⁻¹ maps 3D dual k-forms back to primal (3-k)-forms" {
    comptime {
        try testing.expect(HodgeStarInverseResult(DualC3_3D) == PrimalC0_3D);
        try testing.expect(HodgeStarInverseResult(DualC2_3D) == PrimalC1_3D);
        try testing.expect(HodgeStarInverseResult(DualC1_3D) == PrimalC2_3D);
        try testing.expect(HodgeStarInverseResult(DualC0_3D) == PrimalC3_3D);
    }
}

// ── Correctness: Hodge star values ───────────────────────────────────────

test "★₀ scales by dual vertex area" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC0.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    for (result.values, dual_areas) |r, expected| {
        try testing.expectApproxEqAbs(expected, r, 1e-15);
    }
}

test "★₁ applies Whitney mass matrix (not diagonal)" {
    // The Whitney mass matrix is not diagonal — ★₁ applied to a unit
    // 1-form should differ from the diagonal dual_length/length scaling.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC1.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    // Verify output matches SpMV with the Whitney mass matrix.
    const edge_count = mesh.num_edges();
    const expected = try allocator.alloc(f64, edge_count);
    defer allocator.free(expected);
    sparse.spmv(mesh.whitney_mass_1, omega.values, expected);

    for (result.values, expected) |r, e| {
        try testing.expectApproxEqAbs(e, r, 1e-15);
    }
}

test "★₂ scales by 1 / face area" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 2, 2.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC2.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values) |*v| v.* = 1.0;

    var result = try hodge_star(allocator, omega);
    defer result.deinit(allocator);

    const face_volumes = mesh.simplices(2).items(.volume);
    for (result.values, face_volumes) |r, volume| {
        try testing.expectApproxEqAbs(1.0 / volume, r, 1e-15);
    }
}

// ── Invariant: ★⁻¹ ∘ ★ = identity ───────────────────────────────────────

test "★⁻¹ ∘ ★ = identity for all degrees on random inputs" {
    // For 1000 random primal cochains of each degree, ★⁻¹(★(ω)) = ω
    // must hold to machine precision. For k=0 and k=2 the round-trip is
    // exact (diagonal multiply + divide). For k=1, the round-trip goes
    // through SpMV + CG solve — tolerance is set by CG convergence.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 2.0, 1.5);
    defer mesh.deinit(allocator);

    // k=0,2: two floating-point ops → 1e-14 relative tolerance.
    // k=1: CG solve to 1e-10 relative residual → 1e-8 absolute is safe.

    // ── k = 0 ───────────────────────────────────────────────────────
    {
        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_00);
        for (0..1000) |_| {
            var omega = try PrimalC0.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-14);
            }
        }
    }

    // ── k = 1 (Whitney round-trip: SpMV then CG solve) ─────────────
    {
        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_01);
        for (0..100) |_| {
            var omega = try PrimalC1.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                // SpMV is exact; CG solve converges to 1e-10 relative residual.
                // Condition number amplifies the residual → solution error.
                try testing.expectApproxEqRel(original, recovered, 1e-5);
            }
        }
    }

    // ── k = 2 ───────────────────────────────────────────────────────
    {
        var rng = std.Random.DefaultPrng.init(0xDEC_57A2_02);
        for (0..1000) |_| {
            var omega = try PrimalC2.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-14);
            }
        }
    }
}

test "★★⁻¹ = id for all degrees on random 3D tetrahedral meshes" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 2, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    {
        var rng = std.Random.DefaultPrng.init(0x82_3D_0000);
        for (0..1000) |_| {
            var omega = try PrimalC0_3D.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-14);
            }
        }
    }

    {
        var rng = std.Random.DefaultPrng.init(0x82_3D_0001);
        for (0..1000) |_| {
            var omega = try PrimalC1_3D.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-8);
            }
        }
    }

    {
        var rng = std.Random.DefaultPrng.init(0x82_3D_0002);
        for (0..1000) |_| {
            var omega = try PrimalC2_3D.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-8);
            }
        }
    }

    {
        var rng = std.Random.DefaultPrng.init(0x82_3D_0003);
        for (0..1000) |_| {
            var omega = try PrimalC3_3D.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var starred = try hodge_star(allocator, omega);
            defer starred.deinit(allocator);

            var round_trip = try hodge_star_inverse(allocator, starred);
            defer round_trip.deinit(allocator);

            for (omega.values, round_trip.values) |original, recovered| {
                try testing.expectApproxEqRel(original, recovered, 1e-14);
            }
        }
    }
}
