//! Preconditioned conjugate gradient solver for symmetric positive-definite systems.
//!
//! Solves Ax = b where A is SPD, using the standard PCG algorithm:
//!
//!   r₀ = b − Ax₀
//!   z₀ = M⁻¹r₀          (preconditioner application)
//!   p₀ = z₀
//!   for k = 0, 1, 2, ...
//!     αₖ = ⟨rₖ, zₖ⟩ / ⟨pₖ, Apₖ⟩
//!     xₖ₊₁ = xₖ + αₖ pₖ
//!     rₖ₊₁ = rₖ − αₖ Apₖ
//!     if ‖rₖ₊₁‖ / ‖b‖ < tolerance: converged
//!     zₖ₊₁ = M⁻¹rₖ₊₁
//!     βₖ = ⟨rₖ₊₁, zₖ₊₁⟩ / ⟨rₖ, zₖ⟩
//!     pₖ₊₁ = zₖ₊₁ + βₖ pₖ
//!
//! The preconditioner M⁻¹ is supplied as a function pointer that applies
//! the inverse preconditioner in-place: z = M⁻¹r.

const std = @import("std");
const testing = std.testing;
const sparse = @import("sparse.zig");

/// Result of a CG solve.
pub const SolveResult = struct {
    /// Number of iterations performed.
    iterations: u32,
    /// Final relative residual ‖r‖/‖b‖.
    relative_residual: f64,
    /// Whether the solve converged within the tolerance.
    converged: bool,
};

/// Preconditioner function type.
/// Given input z (which initially holds r), overwrite z with M⁻¹r.
pub const Preconditioner = *const fn (z: []f64, context: *const anyopaque) void;

/// Diagonal preconditioner context — stores the diagonal values.
pub const DiagonalPreconditioner = struct {
    diagonal: []const f64,

    pub fn apply(z: []f64, context: *const anyopaque) void {
        const self: *const DiagonalPreconditioner = @ptrCast(@alignCast(context));
        for (z, self.diagonal) |*zi, di| {
            std.debug.assert(di != 0.0);
            zi.* /= di;
        }
    }
};

/// Solve Ax = b using preconditioned conjugate gradient.
///
/// `x` is both the initial guess and the output. For best performance,
/// initialize x to a reasonable guess (e.g., the solution from the
/// previous timestep). Zero initialization always works.
///
/// `preconditioner` and `preconditioner_context` define M⁻¹. Pass null
/// for both to use no preconditioning (identity preconditioner).
///
/// Returns solve statistics. The solution is written into `x`.
pub fn solve(
    matrix: sparse.CsrMatrix(f64),
    b: []const f64,
    x: []f64,
    tolerance: f64,
    max_iterations: u32,
    preconditioner: ?Preconditioner,
    preconditioner_context: ?*const anyopaque,
    /// Scratch space: 4 vectors of length matrix.n_rows. Caller allocates to avoid
    /// per-solve allocation.
    scratch: Scratch,
) SolveResult {
    const row_count = matrix.n_rows;
    std.debug.assert(b.len == row_count);
    std.debug.assert(x.len == row_count);

    const r = scratch.r;
    const z = scratch.z;
    const p = scratch.p;
    const ap = scratch.ap;

    // r₀ = b − Ax₀
    sparse.spmv(matrix, x, r);
    for (0..row_count) |i| r[i] = b[i] - r[i];

    // Apply preconditioner: z₀ = M⁻¹r₀
    @memcpy(z, r);
    if (preconditioner) |precond| {
        precond(z, preconditioner_context.?);
    }

    // p₀ = z₀
    @memcpy(p, z);

    var rz = dot(r, z);
    const b_norm = @sqrt(dot(b, b));
    if (b_norm == 0.0) return .{ .iterations = 0, .relative_residual = 0.0, .converged = true };

    for (0..max_iterations) |iter| {
        // Check convergence.
        const r_norm = @sqrt(dot(r, r));
        if (r_norm / b_norm < tolerance) {
            return .{
                .iterations = @intCast(iter),
                .relative_residual = r_norm / b_norm,
                .converged = true,
            };
        }

        // αₖ = ⟨rₖ, zₖ⟩ / ⟨pₖ, Apₖ⟩
        sparse.spmv(matrix, p, ap);
        const p_ap = dot(p, ap);
        std.debug.assert(p_ap > 0.0); // A is SPD
        const alpha = rz / p_ap;

        // xₖ₊₁ = xₖ + αpₖ, rₖ₊₁ = rₖ − αApₖ
        for (0..row_count) |i| {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        // zₖ₊₁ = M⁻¹rₖ₊₁
        @memcpy(z, r);
        if (preconditioner) |precond| {
            precond(z, preconditioner_context.?);
        }

        // βₖ = ⟨rₖ₊₁, zₖ₊₁⟩ / ⟨rₖ, zₖ⟩
        const rz_new = dot(r, z);
        const beta = rz_new / rz;
        rz = rz_new;

        // pₖ₊₁ = zₖ₊₁ + βpₖ
        for (0..row_count) |i| p[i] = z[i] + beta * p[i];
    }

    const final_r_norm = @sqrt(dot(r, r));
    return .{
        .iterations = max_iterations,
        .relative_residual = final_r_norm / b_norm,
        .converged = false,
    };
}

/// Scratch space for the CG solver — 4 vectors of length `vector_length`.
pub const Scratch = struct {
    r: []f64,
    z: []f64,
    p: []f64,
    ap: []f64,

    pub fn init(allocator: std.mem.Allocator, vector_length: u32) !Scratch {
        const r = try allocator.alloc(f64, vector_length);
        errdefer allocator.free(r);
        const z = try allocator.alloc(f64, vector_length);
        errdefer allocator.free(z);
        const p = try allocator.alloc(f64, vector_length);
        errdefer allocator.free(p);
        const ap = try allocator.alloc(f64, vector_length);
        return .{ .r = r, .z = z, .p = p, .ap = ap };
    }

    pub fn deinit(self: *Scratch, allocator: std.mem.Allocator) void {
        allocator.free(self.r);
        allocator.free(self.z);
        allocator.free(self.p);
        allocator.free(self.ap);
    }
};

fn dot(a: []const f64, b: []const f64) f64 {
    std.debug.assert(a.len == b.len);
    var sum: f64 = 0;
    for (a, b) |ai, bi| sum += ai * bi;
    return sum;
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

test "CG solves 3×3 diagonal system exactly" {
    const allocator = testing.allocator;
    // A = diag(2, 3, 5), b = [4, 9, 25] → x = [2, 3, 5]
    var m = try sparse.CsrMatrix(f64).init(allocator, 3, 3, 3);
    defer m.deinit(allocator);
    m.row_ptr[0] = 0;
    m.row_ptr[1] = 1;
    m.row_ptr[2] = 2;
    m.row_ptr[3] = 3;
    m.col_idx[0] = 0;
    m.col_idx[1] = 1;
    m.col_idx[2] = 2;
    m.values[0] = 2.0;
    m.values[1] = 3.0;
    m.values[2] = 5.0;

    const b = [_]f64{ 4.0, 9.0, 25.0 };
    var x = [_]f64{ 0.0, 0.0, 0.0 };

    var scratch = try Scratch.init(allocator, 3);
    defer scratch.deinit(allocator);

    const result = solve(m, &b, &x, 1e-12, 100, null, null, scratch);

    try testing.expect(result.converged);
    try testing.expect(result.iterations <= 3);
    try testing.expectApproxEqAbs(@as(f64, 2.0), x[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 3.0), x[1], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 5.0), x[2], 1e-12);
}

test "CG solves dense SPD system" {
    const allocator = testing.allocator;
    // A = [[4, 1], [1, 3]], b = [1, 2] → x = [1/11, 7/11]
    var m = try sparse.CsrMatrix(f64).init(allocator, 2, 2, 4);
    defer m.deinit(allocator);
    m.row_ptr[0] = 0;
    m.row_ptr[1] = 2;
    m.row_ptr[2] = 4;
    m.col_idx[0] = 0;
    m.col_idx[1] = 1;
    m.col_idx[2] = 0;
    m.col_idx[3] = 1;
    m.values[0] = 4.0;
    m.values[1] = 1.0;
    m.values[2] = 1.0;
    m.values[3] = 3.0;

    const b = [_]f64{ 1.0, 2.0 };
    var x = [_]f64{ 0.0, 0.0 };

    var scratch = try Scratch.init(allocator, 2);
    defer scratch.deinit(allocator);

    const result = solve(m, &b, &x, 1e-12, 100, null, null, scratch);

    try testing.expect(result.converged);
    try testing.expectApproxEqAbs(1.0 / 11.0, x[0], 1e-12);
    try testing.expectApproxEqAbs(7.0 / 11.0, x[1], 1e-12);
}

test "CG with diagonal preconditioner converges faster" {
    const allocator = testing.allocator;
    // A = [[100, 1], [1, 1]], b = [101, 2]
    // Poorly conditioned (κ ≈ 100). Diagonal preconditioner helps.
    var m = try sparse.CsrMatrix(f64).init(allocator, 2, 2, 4);
    defer m.deinit(allocator);
    m.row_ptr[0] = 0;
    m.row_ptr[1] = 2;
    m.row_ptr[2] = 4;
    m.col_idx[0] = 0;
    m.col_idx[1] = 1;
    m.col_idx[2] = 0;
    m.col_idx[3] = 1;
    m.values[0] = 100.0;
    m.values[1] = 1.0;
    m.values[2] = 1.0;
    m.values[3] = 1.0;

    const b = [_]f64{ 101.0, 2.0 };
    var x = [_]f64{ 0.0, 0.0 };

    var scratch = try Scratch.init(allocator, 2);
    defer scratch.deinit(allocator);

    var precond = DiagonalPreconditioner{ .diagonal = &[_]f64{ 100.0, 1.0 } };

    const result = solve(m, &b, &x, 1e-12, 100, &precond, scratch);

    try testing.expect(result.converged);
    try testing.expectApproxEqAbs(@as(f64, 1.0), x[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 1.0), x[1], 1e-10);
}

const MockScalePreconditioner = struct {
    scale: f64,

    pub fn apply(self: *const @This(), z: []f64) void {
        std.debug.assert(self.scale != 0.0);
        for (z) |*zi| {
            zi.* /= self.scale;
        }
    }
};

test "PreconditionerConcept accepts a conforming type" {
    comptime PreconditionerConcept(MockScalePreconditioner);
}

test "CG respects the iteration bound when tolerance is not met" {
    const allocator = testing.allocator;

    var matrix = try weighted_path_laplacian(allocator, &[_]f64{ 1.0, 4.0, 9.0, 16.0 });
    defer matrix.deinit(allocator);

    const rhs = [_]f64{ 1.0, 0.0, 0.0, 1.0, 2.0 };
    var x = [_]f64{ 0.0, 0.0, 0.0, 0.0, 0.0 };

    var scratch = try Scratch.init(allocator, matrix.n_rows);
    defer scratch.deinit(allocator);

    const result = solve(matrix, &rhs, &x, 1e-18, 1, null, scratch);

    try testing.expect(!result.converged);
    try testing.expectEqual(@as(u32, 1), result.iterations);
}

test "Jacobi preconditioner reduces iterations on a weighted Laplacian system" {
    const allocator = testing.allocator;

    const edge_weights = [_]f64{ 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0 };
    var matrix = try weighted_path_laplacian(allocator, &edge_weights);
    defer matrix.deinit(allocator);

    const exact_solution = [_]f64{ 0.25, -0.5, 0.75, -1.0, 1.25, -1.5, 1.75 };
    var rhs = [_]f64{0} ** exact_solution.len;
    sparse.spmv(matrix, &exact_solution, &rhs);

    var x_unpreconditioned = [_]f64{0} ** exact_solution.len;
    var x_jacobi = [_]f64{0} ** exact_solution.len;

    var scratch_unpreconditioned = try Scratch.init(allocator, matrix.n_rows);
    defer scratch_unpreconditioned.deinit(allocator);
    var scratch_jacobi = try Scratch.init(allocator, matrix.n_rows);
    defer scratch_jacobi.deinit(allocator);

    const diagonal = try diagonal_of(allocator, matrix);
    defer allocator.free(diagonal);
    var jacobi = DiagonalPreconditioner{ .diagonal = diagonal };

    const unpreconditioned = solve(
        matrix,
        &rhs,
        &x_unpreconditioned,
        1e-10,
        512,
        null,
        scratch_unpreconditioned,
    );
    const preconditioned = solve(
        matrix,
        &rhs,
        &x_jacobi,
        1e-10,
        512,
        &jacobi,
        scratch_jacobi,
    );

    try testing.expect(unpreconditioned.converged);
    try testing.expect(preconditioned.converged);
    try testing.expect(preconditioned.iterations < unpreconditioned.iterations);

    for (x_jacobi, exact_solution) |computed, expected| {
        try testing.expectApproxEqAbs(expected, computed, 1e-8);
    }
}

test "CG solves Whitney mass matrix system" {
    // Integration test: assemble Whitney mass matrix on a small mesh,
    // solve Mx = b for a random b, verify Ax ≈ b.
    const topology = @import("../topology/mesh.zig");
    const whitney = @import("../operators/whitney_mass.zig");

    const allocator = testing.allocator;
    const Mesh2D = topology.Mesh(2, 2);
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try whitney.assemble_whitney_mass(1, allocator, &mesh);
    defer mass.deinit(allocator);

    const row_count = mass.n_rows;

    // Random RHS.
    var rng = std.Random.DefaultPrng.init(0xDEC_C6_00);
    const b_buf = try allocator.alloc(f64, row_count);
    defer allocator.free(b_buf);
    for (b_buf) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    const x_buf = try allocator.alloc(f64, row_count);
    defer allocator.free(x_buf);
    @memset(x_buf, 0.0);

    var scratch = try Scratch.init(allocator, row_count);
    defer scratch.deinit(allocator);

    const result = solve(mass, b_buf, x_buf, 1e-10, 1000, null, scratch);
    try testing.expect(result.converged);

    // Verify: Ax ≈ b.
    const ax = try allocator.alloc(f64, row_count);
    defer allocator.free(ax);
    sparse.spmv(mass, x_buf, ax);

    for (ax, b_buf) |computed, expected| {
        try testing.expectApproxEqAbs(expected, computed, 1e-8);
    }
}

fn weighted_path_laplacian(
    allocator: std.mem.Allocator,
    edge_weights: []const f64,
) !sparse.CsrMatrix(f64) {
    const node_count = edge_weights.len + 1;
    const nonzero_count = 3 * node_count - 2;

    var matrix = try sparse.CsrMatrix(f64).init(
        allocator,
        @intCast(node_count),
        @intCast(node_count),
        @intCast(nonzero_count),
    );
    errdefer matrix.deinit(allocator);

    var diagonal = try allocator.alloc(f64, node_count);
    defer allocator.free(diagonal);
    @memset(diagonal, 0.0);
    for (edge_weights, 0..) |weight, edge_idx| {
        diagonal[edge_idx] += weight;
        diagonal[edge_idx + 1] += weight;
    }

    var nnz_cursor: u32 = 0;
    matrix.row_ptr[0] = 0;
    for (0..node_count) |row_idx| {
        if (row_idx > 0) {
            matrix.col_idx[nnz_cursor] = @intCast(row_idx - 1);
            matrix.values[nnz_cursor] = -edge_weights[row_idx - 1];
            nnz_cursor += 1;
        }

        matrix.col_idx[nnz_cursor] = @intCast(row_idx);
        matrix.values[nnz_cursor] = diagonal[row_idx];
        nnz_cursor += 1;

        if (row_idx + 1 < node_count) {
            matrix.col_idx[nnz_cursor] = @intCast(row_idx + 1);
            matrix.values[nnz_cursor] = -edge_weights[row_idx];
            nnz_cursor += 1;
        }

        matrix.row_ptr[row_idx + 1] = nnz_cursor;
    }

    std.debug.assert(nnz_cursor == nonzero_count);
    return matrix;
}

fn diagonal_of(allocator: std.mem.Allocator, matrix: sparse.CsrMatrix(f64)) ![]f64 {
    const diagonal = try allocator.alloc(f64, matrix.n_rows);
    errdefer allocator.free(diagonal);

    for (0..matrix.n_rows) |row_idx_usize| {
        const row_idx: u32 = @intCast(row_idx_usize);
        diagonal[row_idx] = 0.0;

        const row = matrix.row(row_idx);
        for (row.cols, row.vals) |col_idx, value| {
            if (col_idx != row_idx) continue;
            diagonal[row_idx] = value;
            break;
        }

        std.debug.assert(diagonal[row_idx] > 0.0);
    }

    return diagonal;
}
