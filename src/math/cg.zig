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
    /// Scratch space: 4 vectors of length n. Caller allocates to avoid
    /// per-solve allocation.
    scratch: Scratch,
) SolveResult {
    const n = matrix.n_rows;
    std.debug.assert(b.len == n);
    std.debug.assert(x.len == n);

    const r = scratch.r;
    const z = scratch.z;
    const p = scratch.p;
    const ap = scratch.ap;

    // r₀ = b − Ax₀
    sparse.spmv(matrix, x, r);
    for (0..n) |i| r[i] = b[i] - r[i];

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
        for (0..n) |i| {
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
        for (0..n) |i| p[i] = z[i] + beta * p[i];
    }

    const final_r_norm = @sqrt(dot(r, r));
    return .{
        .iterations = max_iterations,
        .relative_residual = final_r_norm / b_norm,
        .converged = false,
    };
}

/// Scratch space for the CG solver — 4 vectors of length n.
pub const Scratch = struct {
    r: []f64,
    z: []f64,
    p: []f64,
    ap: []f64,

    pub fn init(allocator: std.mem.Allocator, n: u32) !Scratch {
        return .{
            .r = try allocator.alloc(f64, n),
            .z = try allocator.alloc(f64, n),
            .p = try allocator.alloc(f64, n),
            .ap = try allocator.alloc(f64, n),
        };
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

    const result = solve(m, &b, &x, 1e-12, 100, DiagonalPreconditioner.apply, @ptrCast(&precond), scratch);

    try testing.expect(result.converged);
    try testing.expectApproxEqAbs(@as(f64, 1.0), x[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 1.0), x[1], 1e-10);
}

test "CG solves Whitney mass matrix system" {
    // Integration test: assemble Whitney mass matrix on a small mesh,
    // solve Mx = b for a random b, verify Ax ≈ b.
    const topology = @import("../topology/mesh.zig");
    const whitney = @import("../operators/whitney_mass.zig");

    const allocator = testing.allocator;
    const Mesh2D = topology.Mesh(2);
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try whitney.assemble_whitney_mass_1(allocator, &mesh);
    defer mass.deinit(allocator);

    const n = mass.n_rows;

    // Random RHS.
    var rng = std.Random.DefaultPrng.init(0xDEC_C6_00);
    const b_buf = try allocator.alloc(f64, n);
    defer allocator.free(b_buf);
    for (b_buf) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    const x_buf = try allocator.alloc(f64, n);
    defer allocator.free(x_buf);
    @memset(x_buf, 0.0);

    var scratch = try Scratch.init(allocator, n);
    defer scratch.deinit(allocator);

    const result = solve(mass, b_buf, x_buf, 1e-10, 1000, null, null, scratch);
    try testing.expect(result.converged);

    // Verify: Ax ≈ b.
    const ax = try allocator.alloc(f64, n);
    defer allocator.free(ax);
    sparse.spmv(mass, x_buf, ax);

    for (ax, b_buf) |computed, expected| {
        try testing.expectApproxEqAbs(expected, computed, 1e-8);
    }
}
