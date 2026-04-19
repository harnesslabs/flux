//! Discrete Hodge star operator and its inverse.
//!
//! The Hodge star ★ₖ maps primal k-cochains to dual (n−k)-cochains.
//! The implementation dispatches by degree at comptime:
//!
//!   Diagonal path for boundary degrees k = 0 and k = n.
//!   Whitney/Galerkin mass-matrix path for every interior degree 0 < k < n.
//!
//! The mesh stores the interior-degree Whitney mass matrices and diagonal
//! preconditioners, indexed by degree.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const sparse = @import("../math/sparse.zig");
const conjugate_gradient = @import("../math/cg.zig");
const whitney_mass = @import("whitney_mass.zig");

pub const SolveError = error{
    HodgeStarInverseDidNotConverge,
};

pub const MetricMode = enum {
    flat,
    riemannian,
};

pub const MetricError = error{
    MetricNotYetImplemented,
    UnsupportedMetricDegree,
};

pub fn Metric(comptime MeshType: type, comptime mode: MetricMode) type {
    return switch (mode) {
        .flat => struct {
            pub const metric_mode = MetricMode.flat;
        },
        .riemannian => struct {
            pub const metric_mode = MetricMode.riemannian;

            top_simplex_tensors: []const [MeshType.topological_dimension][MeshType.topological_dimension]f64,
        },
    };
}

const WhitneyInverseSolveParams = struct {
    relative_tolerance: f64,
    iteration_limit: u32,
};

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
/// For k=0 and k=n, this is a diagonal scaling by dual/primal measures.
/// For 0 < k < n, this applies the Whitney/Galerkin mass matrix Mₖ via SpMV.
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

    if (comptime supportsWhitneyMassDegree(InputType.MeshT, k)) {
        sparse.spmv(input.mesh.whitney_mass(k), input.values, output.values);
        return output;
    }

    applyDiagonalForward(InputType.MeshT, k, input.mesh, input.values, output.values);
    return output;
}

/// Apply the inverse Hodge star ★⁻¹ to a dual (n−k)-cochain, returning a
/// primal k-cochain.
///
/// For k=0 and k=n, this is the element-wise reciprocal of the diagonal.
/// For 0 < k < n, this solves Mₖ x = b via preconditioned conjugate gradient.
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

    if (comptime supportsWhitneyMassDegree(InputType.MeshT, primal_degree)) {
        const solve_params = whitneyInverseSolveParams(InputType.MeshT, primal_degree);
        try solveWhitneyInverse(
            allocator,
            input.mesh.whitney_mass(primal_degree),
            input.mesh.whitney_preconditioner(primal_degree),
            input.values,
            output.values,
            solve_params.relative_tolerance,
            solve_params.iteration_limit,
        );
        return output;
    }

    try applyDiagonalInverse(allocator, InputType.MeshT, primal_degree, input.mesh, input.values, output.values);
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

pub fn apply_hodge_star_with_metric(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    return switch (@TypeOf(metric).metric_mode) {
        .flat => apply_hodge_star(allocator, input),
        .riemannian => apply_hodge_star_riemannian(allocator, metric, input),
    };
}

pub fn apply_hodge_star_inverse_with_metric(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    return switch (@TypeOf(metric).metric_mode) {
        .flat => apply_hodge_star_inverse(allocator, input),
        .riemannian => apply_hodge_star_inverse_riemannian(allocator, metric, input),
    };
}

pub fn hodge_star_with_metric(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    return apply_hodge_star_with_metric(allocator, metric, input);
}

pub fn hodge_star_inverse_with_metric(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    return apply_hodge_star_inverse_with_metric(allocator, metric, input);
}

fn apply_hodge_star_riemannian(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInput(InputType);

    const MeshType = InputType.MeshT;
    const k = InputType.degree;
    const n = MeshType.topological_dimension;
    std.debug.assert(metric.top_simplex_tensors.len == input.mesh.num_cells(n));

    var output = try HodgeStarResult(InputType).init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    if (comptime supportsWhitneyMassDegree(MeshType, k)) {
        var matrix = try whitney_mass.assemble_whitney_mass(k, allocator, input.mesh, metric.top_simplex_tensors);
        defer matrix.deinit(allocator);
        sparse.spmv(matrix, input.values, output.values);
        return output;
    }

    if (k == n) {
        const primal_volumes = input.mesh.simplices(n).items(.volume);
        for (output.values, input.values, primal_volumes, metric.top_simplex_tensors) |*out, in_value, volume, tensor| {
            const metric_volume = volume * @sqrt(metricTensorDeterminant(n, tensor));
            std.debug.assert(metric_volume != 0.0);
            out.* = in_value / metric_volume;
        }
        return output;
    }

    return MetricError.UnsupportedMetricDegree;
}

fn apply_hodge_star_inverse_riemannian(
    allocator: std.mem.Allocator,
    metric: anytype,
    input: anytype,
) !HodgeStarInverseResult(@TypeOf(input)) {
    const InputType = @TypeOf(input);
    comptime validateHodgeStarInverseInput(InputType);

    const MeshType = InputType.MeshT;
    const n = MeshType.topological_dimension;
    const primal_degree = n - InputType.degree;
    std.debug.assert(metric.top_simplex_tensors.len == input.mesh.num_cells(n));

    var output = try HodgeStarInverseResult(InputType).init(allocator, input.mesh);
    errdefer output.deinit(allocator);

    if (comptime supportsWhitneyMassDegree(MeshType, primal_degree)) {
        var matrix = try whitney_mass.assemble_whitney_mass(primal_degree, allocator, input.mesh, metric.top_simplex_tensors);
        defer matrix.deinit(allocator);

        const diagonal = try assembleMatrixDiagonal(allocator, matrix);
        defer allocator.free(diagonal);

        try solveWhitneyInverse(
            allocator,
            matrix,
            diagonal,
            input.values,
            output.values,
            whitneyInverseSolveParams(MeshType, primal_degree).relative_tolerance,
            whitneyInverseSolveParams(MeshType, primal_degree).iteration_limit,
        );
        return output;
    }

    if (primal_degree == n) {
        const primal_volumes = input.mesh.simplices(n).items(.volume);
        for (output.values, input.values, primal_volumes, metric.top_simplex_tensors) |*out, in_value, volume, tensor| {
            out.* = in_value * volume * @sqrt(metricTensorDeterminant(n, tensor));
        }
        return output;
    }

    return MetricError.UnsupportedMetricDegree;
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
    if (comptime supportsWhitneyMassDegree(MeshType, primal_degree)) {
        sparse.spmv(mesh.whitney_mass(primal_degree), input, output);
        return;
    }

    applyDiagonalForward(MeshType, primal_degree, mesh, input, output);
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
    if (comptime supportsWhitneyMassDegree(MeshType, primal_degree)) {
        const solve_params = whitneyInverseSolveParams(MeshType, primal_degree);
        try solveWhitneyInverse(
            allocator,
            mesh.whitney_mass(primal_degree),
            mesh.whitney_preconditioner(primal_degree),
            input,
            output,
            solve_params.relative_tolerance,
            solve_params.iteration_limit,
        );
        return;
    }

    try applyDiagonalInverse(allocator, MeshType, primal_degree, mesh, input, output);
}

fn solveWhitneyInverse(
    allocator: std.mem.Allocator,
    matrix: sparse.CsrMatrix(f64),
    diagonal: []const f64,
    input: []const f64,
    output: []f64,
    relative_tolerance: f64,
    iteration_limit: u32,
) !void {
    std.debug.assert(diagonal.len == output.len);

    // Use the diagonal DEC ratio as both preconditioner and initial guess.
    for (output, input, diagonal) |*out, in_val, diagonal_entry| {
        std.debug.assert(diagonal_entry != 0.0);
        out.* = in_val / diagonal_entry;
    }

    var scratch = try conjugate_gradient.Scratch.init(allocator, @intCast(output.len));
    defer scratch.deinit(allocator);

    var precond = conjugate_gradient.DiagonalPreconditioner{ .diagonal = diagonal };
    const result = conjugate_gradient.solve(
        matrix,
        input,
        output,
        relative_tolerance,
        iteration_limit,
        &precond,
        scratch,
    );
    if (!result.converged) {
        return SolveError.HodgeStarInverseDidNotConverge;
    }
}

fn assembleMatrixDiagonal(
    allocator: std.mem.Allocator,
    matrix: sparse.CsrMatrix(f64),
) ![]f64 {
    const diagonal = try allocator.alloc(f64, matrix.n_rows);
    errdefer allocator.free(diagonal);
    @memset(diagonal, 0.0);

    for (0..matrix.n_rows) |row_idx| {
        const row = matrix.row(@intCast(row_idx));
        for (row.cols, row.vals) |col_idx, value| {
            if (col_idx != row_idx) continue;
            diagonal[row_idx] += value;
        }
        std.debug.assert(diagonal[row_idx] > 0.0);
    }

    return diagonal;
}

// ── Diagonal application ────────────────────────────────────────────────

fn applyDiagonalForward(
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
) void {
    std.debug.assert(!supportsWhitneyMassDegree(MeshType, primal_degree));

    if (primal_degree == 0) {
        const dual_volumes = mesh.vertices.slice().items(.dual_volume);
        for (output, input, dual_volumes) |*out, in_val, ratio| {
            out.* = ratio * in_val;
        }
        return;
    }

    if (primal_degree == MeshType.topological_dimension) {
        const primal_volumes = mesh.simplices(MeshType.topological_dimension).items(.volume);
        for (output, input, primal_volumes) |*out, in_val, volume| {
            std.debug.assert(volume != 0.0);
            out.* = in_val / volume;
        }
        return;
    }

    unreachable;
}

fn applyDiagonalInverse(
    allocator: std.mem.Allocator,
    comptime MeshType: type,
    comptime primal_degree: comptime_int,
    mesh: *const MeshType,
    input: []const f64,
    output: []f64,
) !void {
    std.debug.assert(!supportsWhitneyMassDegree(MeshType, primal_degree));
    _ = allocator;

    if (primal_degree == 0) {
        const dual_volumes = mesh.vertices.slice().items(.dual_volume);
        for (output, input, dual_volumes) |*out, in_val, ratio| {
            std.debug.assert(ratio != 0.0);
            out.* = in_val / ratio;
        }
        return;
    }

    if (primal_degree == MeshType.topological_dimension) {
        const primal_volumes = mesh.simplices(MeshType.topological_dimension).items(.volume);
        for (output, input, primal_volumes) |*out, in_val, volume| {
            out.* = volume * in_val;
        }
        return;
    }

    unreachable;
}

fn supportsWhitneyMassDegree(comptime MeshType: type, comptime primal_degree: comptime_int) bool {
    return primal_degree > 0 and primal_degree < MeshType.topological_dimension;
}

fn metricTensorDeterminant(
    comptime n: comptime_int,
    tensor: [n][n]f64,
) f64 {
    if (n == 2) {
        return tensor[0][0] * tensor[1][1] - tensor[0][1] * tensor[1][0];
    }

    if (n == 3) {
        const a = tensor[0][0];
        const b = tensor[0][1];
        const c = tensor[0][2];
        const d = tensor[1][0];
        const e = tensor[1][1];
        const f = tensor[1][2];
        const g = tensor[2][0];
        const h = tensor[2][1];
        const i = tensor[2][2];
        return a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g);
    }

    @compileError("metric tensor determinant is only implemented for n = 2 or n = 3");
}

fn whitneyInverseSolveParams(comptime MeshType: type, comptime primal_degree: comptime_int) WhitneyInverseSolveParams {
    std.debug.assert(supportsWhitneyMassDegree(MeshType, primal_degree));

    if (MeshType.topological_dimension <= 2) {
        // Preserve the original 2D benchmarked solve target so the PR remains
        // comparable to main on the existing hodge_star_inverse_1_cg row.
        return .{
            .relative_tolerance = 1e-10,
            .iteration_limit = 1000,
        };
    }

    // Random 3D tetrahedral meshes need a tighter residual target for the
    // ★★⁻¹ identity property to hold across distorted cells.
    return .{
        .relative_tolerance = 1e-14,
        .iteration_limit = 20000,
    };
}

fn expectApproxEqRelOrAbs(expected: f64, actual: f64, relative_tolerance: f64, absolute_tolerance: f64) !void {
    const tolerance = @max(absolute_tolerance, @abs(expected) * relative_tolerance);
    try testing.expect(@abs(expected - actual) <= tolerance);
}

fn expectSlicesApproxEqAbs(expected: []const f64, actual: []const f64, tolerance: f64) !void {
    try testing.expectEqual(expected.len, actual.len);
    for (expected, actual) |expected_value, actual_value| {
        try testing.expectApproxEqAbs(expected_value, actual_value, tolerance);
    }
}

fn expectRoundTripIdentity3DForDegree(
    comptime degree: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const Mesh3D,
    random: std.Random,
    trial_count: u32,
    relative_tolerance: f64,
    absolute_tolerance: f64,
) !void {
    const PrimalType = cochain.Cochain(Mesh3D, degree, cochain.Primal);

    for (0..trial_count) |_| {
        var omega = try PrimalType.init(allocator, mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*value| value.* = random.float(f64) * 200.0 - 100.0;

        var starred = try hodge_star(allocator, omega);
        defer starred.deinit(allocator);

        var round_trip = try hodge_star_inverse(allocator, starred);
        defer round_trip.deinit(allocator);

        for (omega.values, round_trip.values) |original, recovered| {
            try expectApproxEqRelOrAbs(original, recovered, relative_tolerance, absolute_tolerance);
        }
    }
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
const MeshSurface = topology.Mesh(3, 2);

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
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
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
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
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
    sparse.spmv(mesh.whitney_mass(1), omega.values, expected);

    for (result.values, expected) |r, e| {
        try testing.expectApproxEqAbs(e, r, 1e-15);
    }
}

test "★₂ scales by 1 / face area" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 2, 2.0, 1.0);
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

test "Metric(.flat) reproduces Euclidean Hodge star exactly for all 2D degrees" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 2, 2.0, 1.5);
    defer mesh.deinit(allocator);

    const metric = Metric(Mesh2D, .flat){};

    {
        var omega = try PrimalC0.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values, 0..) |*value, idx| value.* = @as(f64, @floatFromInt(idx + 1)) * 0.5;

        var baseline = try hodge_star(allocator, omega);
        defer baseline.deinit(allocator);

        var metric_result = try hodge_star_with_metric(allocator, metric, omega);
        defer metric_result.deinit(allocator);

        try expectSlicesApproxEqAbs(baseline.values, metric_result.values, 1e-15);
    }

    {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values, 0..) |*value, idx| value.* = @as(f64, @floatFromInt(idx + 1)) * -0.25;

        var baseline = try hodge_star(allocator, omega);
        defer baseline.deinit(allocator);

        var metric_result = try hodge_star_with_metric(allocator, metric, omega);
        defer metric_result.deinit(allocator);

        try expectSlicesApproxEqAbs(baseline.values, metric_result.values, 1e-15);
    }

    {
        var omega = try PrimalC2.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values, 0..) |*value, idx| value.* = @as(f64, @floatFromInt(idx + 1)) * 0.75;

        var baseline = try hodge_star(allocator, omega);
        defer baseline.deinit(allocator);

        var metric_result = try hodge_star_with_metric(allocator, metric, omega);
        defer metric_result.deinit(allocator);

        try expectSlicesApproxEqAbs(baseline.values, metric_result.values, 1e-15);
    }
}

test "Metric(.riemannian) on a constant tensor matches Euclidean star on the pulled-back mesh" {
    const allocator = testing.allocator;
    const faces = [_][3]u32{
        .{ 0, 1, 2 },
        .{ 0, 2, 3 },
    };
    const euclidean_vertices = [_][2]f64{
        .{ 0.0, 0.0 },
        .{ 1.0, 0.0 },
        .{ 1.0, 1.0 },
        .{ 0.0, 1.0 },
    };
    const transformed_vertices = [_][2]f64{
        .{ 0.0, 0.0 },
        .{ 2.0, 0.0 },
        .{ 2.0, 1.0 },
        .{ 0.0, 1.0 },
    };

    var base_mesh = try Mesh2D.from_triangles(allocator, &euclidean_vertices, &faces);
    defer base_mesh.deinit(allocator);

    var transformed_mesh = try Mesh2D.from_triangles(allocator, &transformed_vertices, &faces);
    defer transformed_mesh.deinit(allocator);

    const tensors = [_][2][2]f64{
        .{ .{ 4.0, 0.0 }, .{ 0.0, 1.0 } },
        .{ .{ 4.0, 0.0 }, .{ 0.0, 1.0 } },
    };
    const metric = Metric(Mesh2D, .riemannian){
        .top_simplex_tensors = &tensors,
    };

    var omega_base = try PrimalC1.init(allocator, &base_mesh);
    defer omega_base.deinit(allocator);
    var omega_transformed = try PrimalC1.init(allocator, &transformed_mesh);
    defer omega_transformed.deinit(allocator);

    for (omega_base.values, 0..) |*value, idx| {
        value.* = @as(f64, @floatFromInt(idx + 1)) * 0.25;
        omega_transformed.values[idx] = value.*;
    }

    var metric_result = try hodge_star_with_metric(allocator, metric, omega_base);
    defer metric_result.deinit(allocator);

    var transformed_result = try hodge_star(allocator, omega_transformed);
    defer transformed_result.deinit(allocator);

    try expectSlicesApproxEqAbs(transformed_result.values, metric_result.values, 1e-12);
}

test "Metric(.riemannian) preserves ★⁻¹ ∘ ★ = id for 2D primal 1-forms" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const face_count = mesh.num_faces();
    const tensors = try allocator.alloc([2][2]f64, face_count);
    defer allocator.free(tensors);
    for (tensors, 0..) |*tensor, face_idx| {
        if ((face_idx % 2) == 0) {
            tensor.* = .{ .{ 2.0, 0.25 }, .{ 0.25, 1.5 } };
        } else {
            tensor.* = .{ .{ 1.5, -0.1 }, .{ -0.1, 1.25 } };
        }
    }

    const metric = Metric(Mesh2D, .riemannian){
        .top_simplex_tensors = tensors,
    };

    var rng = std.Random.DefaultPrng.init(0x85_600D);
    for (0..100) |_| {
        var omega = try PrimalC1.init(allocator, &mesh);
        defer omega.deinit(allocator);
        for (omega.values) |*value| value.* = rng.random().float(f64) * 20.0 - 10.0;

        var starred = try hodge_star_with_metric(allocator, metric, omega);
        defer starred.deinit(allocator);

        var round_trip = try hodge_star_inverse_with_metric(allocator, metric, starred);
        defer round_trip.deinit(allocator);

        for (omega.values, round_trip.values) |expected, actual| {
            try testing.expectApproxEqRel(expected, actual, 1e-6);
        }
    }
}

test "Mesh(3, 2) induced metric matches intrinsic 2D Hodge star on an isometric embedding" {
    const allocator = testing.allocator;
    const faces = [_][3]u32{
        .{ 0, 1, 2 },
        .{ 0, 2, 3 },
    };
    const intrinsic_vertices = [_][2]f64{
        .{ 0.0, 0.0 },
        .{ 1.0, 0.0 },
        .{ 2.0, 1.0 },
        .{ 1.0, 1.0 },
    };
    const embedded_vertices = [_][3]f64{
        .{ 0.0, 0.0, 0.0 },
        .{ 0.7071067811865475, 0.7071067811865475, 0.0 },
        .{ 1.414213562373095, 1.414213562373095, 1.0 },
        .{ 0.7071067811865475, 0.7071067811865475, 1.0 },
    };

    var intrinsic_mesh = try Mesh2D.from_triangles(allocator, &intrinsic_vertices, &faces);
    defer intrinsic_mesh.deinit(allocator);

    var embedded_mesh = try MeshSurface.from_triangles(allocator, &embedded_vertices, &faces);
    defer embedded_mesh.deinit(allocator);

    {
        const EmbeddedPrimalC0 = cochain.Cochain(MeshSurface, 0, cochain.Primal);
        var intrinsic_form = try PrimalC0.init(allocator, &intrinsic_mesh);
        defer intrinsic_form.deinit(allocator);
        var embedded_form = try EmbeddedPrimalC0.init(allocator, &embedded_mesh);
        defer embedded_form.deinit(allocator);

        for (intrinsic_form.values, 0..) |*value, idx| {
            value.* = @as(f64, @floatFromInt(idx + 1)) * 0.5;
            embedded_form.values[idx] = value.*;
        }

        var intrinsic_star = try hodge_star(allocator, intrinsic_form);
        defer intrinsic_star.deinit(allocator);
        var embedded_star = try hodge_star(allocator, embedded_form);
        defer embedded_star.deinit(allocator);

        try expectSlicesApproxEqAbs(intrinsic_star.values, embedded_star.values, 1e-12);
    }

    {
        const EmbeddedPrimalC1 = cochain.Cochain(MeshSurface, 1, cochain.Primal);
        var intrinsic_form = try PrimalC1.init(allocator, &intrinsic_mesh);
        defer intrinsic_form.deinit(allocator);
        var embedded_form = try EmbeddedPrimalC1.init(allocator, &embedded_mesh);
        defer embedded_form.deinit(allocator);

        for (intrinsic_form.values, 0..) |*value, idx| {
            value.* = @as(f64, @floatFromInt(idx + 2)) * -0.25;
            embedded_form.values[idx] = value.*;
        }

        var intrinsic_star = try hodge_star(allocator, intrinsic_form);
        defer intrinsic_star.deinit(allocator);
        var embedded_star = try hodge_star(allocator, embedded_form);
        defer embedded_star.deinit(allocator);

        try expectSlicesApproxEqAbs(intrinsic_star.values, embedded_star.values, 1e-12);
    }

    {
        const EmbeddedPrimalC2 = cochain.Cochain(MeshSurface, 2, cochain.Primal);
        var intrinsic_form = try PrimalC2.init(allocator, &intrinsic_mesh);
        defer intrinsic_form.deinit(allocator);
        var embedded_form = try EmbeddedPrimalC2.init(allocator, &embedded_mesh);
        defer embedded_form.deinit(allocator);

        for (intrinsic_form.values, 0..) |*value, idx| {
            value.* = @as(f64, @floatFromInt(idx + 3)) * 0.75;
            embedded_form.values[idx] = value.*;
        }

        var intrinsic_star = try hodge_star(allocator, intrinsic_form);
        defer intrinsic_star.deinit(allocator);
        var embedded_star = try hodge_star(allocator, embedded_form);
        defer embedded_star.deinit(allocator);

        try expectSlicesApproxEqAbs(intrinsic_star.values, embedded_star.values, 1e-12);
    }
}

// ── Invariant: ★⁻¹ ∘ ★ = identity ───────────────────────────────────────

test "★⁻¹ ∘ ★ = identity for all degrees on random inputs" {
    // For 1000 random primal cochains of each degree, ★⁻¹(★(ω)) = ω
    // must hold to machine precision. For k=0 and k=2 the round-trip is
    // exact (diagonal multiply + divide). For k=1, the round-trip goes
    // through SpMV + CG solve — tolerance is set by CG convergence.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 5, 4, 2.0, 1.5);
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
    // Run 1000 random 3D trials total. Spread them across multiple
    // tetrahedral grids so the property test exercises both coefficient
    // recovery and mesh-dependent orientation/geometry paths.
    const mesh_count: u32 = 10;
    const trials_per_mesh_per_degree: u32 = 25;

    var mesh_rng = std.Random.DefaultPrng.init(0x823D_AE50);
    for (0..mesh_count) |_| {
        const nx: u32 = @intCast(mesh_rng.random().intRangeAtMost(u32, 1, 3));
        const ny: u32 = @intCast(mesh_rng.random().intRangeAtMost(u32, 1, 3));
        const nz: u32 = @intCast(mesh_rng.random().intRangeAtMost(u32, 1, 3));
        const width = 0.5 + mesh_rng.random().float(f64) * 1.5;
        const height = 0.5 + mesh_rng.random().float(f64) * 1.5;
        const depth = 0.5 + mesh_rng.random().float(f64) * 1.5;

        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, nx, ny, nz, width, height, depth);
        defer mesh.deinit(allocator);

        var rng_0 = std.Random.DefaultPrng.init(mesh_rng.random().int(u64));
        try expectRoundTripIdentity3DForDegree(0, allocator, &mesh, rng_0.random(), trials_per_mesh_per_degree, 1e-14, 0.0);

        var rng_1 = std.Random.DefaultPrng.init(mesh_rng.random().int(u64));
        try expectRoundTripIdentity3DForDegree(1, allocator, &mesh, rng_1.random(), trials_per_mesh_per_degree, 1e-7, 1e-10);

        var rng_2 = std.Random.DefaultPrng.init(mesh_rng.random().int(u64));
        try expectRoundTripIdentity3DForDegree(2, allocator, &mesh, rng_2.random(), trials_per_mesh_per_degree, 1e-7, 1e-10);

        var rng_3 = std.Random.DefaultPrng.init(mesh_rng.random().int(u64));
        try expectRoundTripIdentity3DForDegree(3, allocator, &mesh, rng_3.random(), trials_per_mesh_per_degree, 1e-14, 0.0);
    }
}

test "★⁻¹ returns error when Whitney CG solve exhausts iteration limit" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var omega = try PrimalC1.init(allocator, &mesh);
    defer omega.deinit(allocator);
    for (omega.values, 0..) |*value, idx| {
        value.* = @as(f64, @floatFromInt(idx + 1));
    }

    var starred = try hodge_star(allocator, omega);
    defer starred.deinit(allocator);

    var output = try PrimalC1.init(allocator, &mesh);
    defer output.deinit(allocator);

    try testing.expectError(
        SolveError.HodgeStarInverseDidNotConverge,
        solveWhitneyInverse(
            allocator,
            mesh.whitney_mass(1),
            mesh.whitney_preconditioner(1),
            starred.values,
            output.values,
            whitneyInverseSolveParams(Mesh2D, 1).relative_tolerance,
            0,
        ),
    );
}
