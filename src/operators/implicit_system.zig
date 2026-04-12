//! Reusable assembled implicit linear solves with owned workspace.

const std = @import("std");
const testing = std.testing;
const sparse = @import("../math/sparse.zig");
const cg = @import("../math/cg.zig");

pub const SolveConfig = struct {
    tolerance_relative: f64 = 1e-10,
    iteration_limit: u32 = 4000,
};

pub const AssembledImplicitSystem = struct {
    matrix: sparse.CsrMatrix(f64),
    diagonal: []f64,
    rhs: []f64,
    solution: []f64,
    scratch: cg.Scratch,
    config: SolveConfig,

    pub fn init(
        allocator: std.mem.Allocator,
        matrix: sparse.CsrMatrix(f64),
        config: SolveConfig,
    ) !AssembledImplicitSystem {
        std.debug.assert(matrix.n_rows == matrix.n_cols);

        const diagonal = try diagonalOf(allocator, matrix);
        errdefer allocator.free(diagonal);

        const rhs = try allocator.alloc(f64, matrix.n_rows);
        errdefer allocator.free(rhs);
        @memset(rhs, 0.0);

        const solution = try allocator.alloc(f64, matrix.n_rows);
        errdefer allocator.free(solution);
        @memset(solution, 0.0);

        var scratch = try cg.Scratch.init(allocator, @intCast(matrix.n_rows));
        errdefer scratch.deinit(allocator);

        return .{
            .matrix = matrix,
            .diagonal = diagonal,
            .rhs = rhs,
            .solution = solution,
            .scratch = scratch,
            .config = config,
        };
    }

    pub fn deinit(self: *AssembledImplicitSystem, allocator: std.mem.Allocator) void {
        self.scratch.deinit(allocator);
        allocator.free(self.solution);
        allocator.free(self.rhs);
        allocator.free(self.diagonal);
        self.matrix.deinit(allocator);
    }

    pub fn rhsValues(self: *AssembledImplicitSystem) []f64 {
        return self.rhs;
    }

    pub fn solutionValues(self: *AssembledImplicitSystem) []f64 {
        return self.solution;
    }

    pub fn seedSolution(self: *AssembledImplicitSystem, values: []const f64) void {
        std.debug.assert(values.len == self.solution.len);
        @memcpy(self.solution, values);
    }

    pub fn solve(self: *AssembledImplicitSystem) !cg.SolveResult {
        var preconditioner = cg.DiagonalPreconditioner{ .diagonal = self.diagonal };
        const result = cg.solve(
            self.matrix,
            self.rhs,
            self.solution,
            self.config.tolerance_relative,
            self.config.iteration_limit,
            &preconditioner,
            self.scratch,
        );
        if (!result.converged) return error.ConjugateGradientDidNotConverge;
        return result;
    }
};

pub const DirichletConstrainedSystem = struct {
    pub fn init(
        allocator: std.mem.Allocator,
        full_matrix: sparse.CsrMatrix(f64),
        constrained_mask: []const bool,
        config: SolveConfig,
    ) !DirichletConstrainedSystem {
        _ = allocator;
        _ = full_matrix;
        _ = constrained_mask;
        _ = config;
        @panic("not yet implemented");
    }

    pub fn deinit(self: *DirichletConstrainedSystem, allocator: std.mem.Allocator) void {
        _ = self;
        _ = allocator;
        @panic("not yet implemented");
    }

    pub fn fullRhsValues(self: *DirichletConstrainedSystem) []f64 {
        _ = self;
        @panic("not yet implemented");
    }

    pub fn seedSolutionFromFull(self: *DirichletConstrainedSystem, full_values: []const f64) void {
        _ = self;
        _ = full_values;
        @panic("not yet implemented");
    }

    pub fn solveDirichlet(
        self: *DirichletConstrainedSystem,
        boundary_values: []const f64,
        full_solution: []f64,
    ) !cg.SolveResult {
        _ = self;
        _ = boundary_values;
        _ = full_solution;
        @panic("not yet implemented");
    }

    pub fn solveHomogeneousDirichlet(
        self: *DirichletConstrainedSystem,
        full_solution: []f64,
    ) !cg.SolveResult {
        _ = self;
        _ = full_solution;
        @panic("not yet implemented");
    }
};

fn diagonalOf(allocator: std.mem.Allocator, matrix: sparse.CsrMatrix(f64)) ![]f64 {
    const diagonal = try allocator.alloc(f64, matrix.n_rows);
    errdefer allocator.free(diagonal);
    @memset(diagonal, 0.0);

    for (0..matrix.n_rows) |row_idx_usize| {
        const row_idx: u32 = @intCast(row_idx_usize);
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

test "assembled implicit system computes matrix diagonal and reuses owned buffers across solves" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 2, 2, 2);
    matrix.row_ptr[0] = 0;
    matrix.row_ptr[1] = 1;
    matrix.row_ptr[2] = 2;
    matrix.col_idx[0] = 0;
    matrix.col_idx[1] = 1;
    matrix.values[0] = 2.0;
    matrix.values[1] = 4.0;

    var system = try AssembledImplicitSystem.init(allocator, matrix, .{});
    defer system.deinit(allocator);

    const rhs_ptr = system.rhsValues().ptr;
    const solution_ptr = system.solutionValues().ptr;

    system.rhsValues()[0] = 4.0;
    system.rhsValues()[1] = 8.0;
    @memset(system.solutionValues(), 0.0);

    const first_result = try system.solve();
    try testing.expect(first_result.converged);
    try testing.expectApproxEqAbs(@as(f64, 2.0), system.solutionValues()[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 2.0), system.solutionValues()[1], 1e-12);

    try testing.expect(rhs_ptr == system.rhsValues().ptr);
    try testing.expect(solution_ptr == system.solutionValues().ptr);

    system.rhsValues()[0] = 2.0;
    system.rhsValues()[1] = 4.0;

    const second_result = try system.solve();
    try testing.expect(second_result.converged);
    try testing.expectApproxEqAbs(@as(f64, 1.0), system.solutionValues()[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 1.0), system.solutionValues()[1], 1e-12);
}

test "Dirichlet constrained system scatters boundary values and corrects reduced rhs" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 3, 3, 7);
    matrix.row_ptr[0] = 0;
    matrix.row_ptr[1] = 2;
    matrix.row_ptr[2] = 5;
    matrix.row_ptr[3] = 7;
    matrix.col_idx[0] = 0;
    matrix.col_idx[1] = 1;
    matrix.col_idx[2] = 0;
    matrix.col_idx[3] = 1;
    matrix.col_idx[4] = 2;
    matrix.col_idx[5] = 1;
    matrix.col_idx[6] = 2;
    matrix.values[0] = 2.0;
    matrix.values[1] = -1.0;
    matrix.values[2] = -1.0;
    matrix.values[3] = 2.0;
    matrix.values[4] = -1.0;
    matrix.values[5] = -1.0;
    matrix.values[6] = 2.0;

    var system = try DirichletConstrainedSystem.init(
        allocator,
        matrix,
        &[_]bool{ true, false, true },
        .{},
    );
    defer system.deinit(allocator);

    system.fullRhsValues()[0] = 0.0;
    system.fullRhsValues()[1] = 0.0;
    system.fullRhsValues()[2] = 0.0;

    const boundary_values = [_]f64{ 10.0, 0.0, 20.0 };
    var solution = [_]f64{ 0.0, 0.0, 0.0 };

    const result = try system.solveDirichlet(&boundary_values, &solution);
    try testing.expect(result.converged);
    try testing.expectApproxEqAbs(@as(f64, 10.0), solution[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 15.0), solution[1], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 20.0), solution[2], 1e-12);
}

test "Dirichlet constrained system reuses full rhs buffer and seeds from full state" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 3, 3, 3);
    matrix.row_ptr[0] = 0;
    matrix.row_ptr[1] = 1;
    matrix.row_ptr[2] = 2;
    matrix.row_ptr[3] = 3;
    matrix.col_idx[0] = 0;
    matrix.col_idx[1] = 1;
    matrix.col_idx[2] = 2;
    matrix.values[0] = 4.0;
    matrix.values[1] = 5.0;
    matrix.values[2] = 6.0;

    var system = try DirichletConstrainedSystem.init(
        allocator,
        matrix,
        &[_]bool{ true, false, false },
        .{},
    );
    defer system.deinit(allocator);

    const rhs_ptr = system.fullRhsValues().ptr;
    system.seedSolutionFromFull(&[_]f64{ 100.0, 10.0, 18.0 });
    system.fullRhsValues()[0] = 0.0;
    system.fullRhsValues()[1] = 10.0;
    system.fullRhsValues()[2] = 18.0;

    var solution = [_]f64{ 0.0, 0.0, 0.0 };
    const result = try system.solveHomogeneousDirichlet(&solution);
    try testing.expect(result.converged);
    try testing.expect(rhs_ptr == system.fullRhsValues().ptr);
    try testing.expectApproxEqAbs(@as(f64, 0.0), solution[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 2.0), solution[1], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 3.0), solution[2], 1e-12);
}

test {
    testing.refAllDeclsRecursive(@This());
}
