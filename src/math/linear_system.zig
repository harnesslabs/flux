//! Reusable linear-system runtimes and elimination maps.

const std = @import("std");
const testing = std.testing;
const sparse = @import("sparse.zig");
const cg = @import("cg.zig");

pub const SolveConfig = struct {
    tolerance_relative: f64 = 1e-10,
    iteration_limit: u32 = 4000,
};

pub const LinearSystem = struct {
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
    ) !LinearSystem {
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

    pub fn deinit(self: *LinearSystem, allocator: std.mem.Allocator) void {
        self.scratch.deinit(allocator);
        allocator.free(self.solution);
        allocator.free(self.rhs);
        allocator.free(self.diagonal);
        self.matrix.deinit(allocator);
    }

    pub fn rhsValues(self: *LinearSystem) []f64 {
        return self.rhs;
    }

    pub fn solutionValues(self: *LinearSystem) []f64 {
        return self.solution;
    }

    pub fn seedSolution(self: *LinearSystem, values: []const f64) void {
        std.debug.assert(values.len == self.solution.len);
        @memcpy(self.solution, values);
    }

    pub fn solve(self: *LinearSystem) !cg.SolveResult {
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

pub const EliminationMap = struct {
    constrained_mask: []bool,
    reduced_index: []u32,
    free_indices: []u32,

    pub fn init(allocator: std.mem.Allocator, constrained_mask: []const bool) !EliminationMap {
        const owned_mask = try allocator.dupe(bool, constrained_mask);
        errdefer allocator.free(owned_mask);

        const reduced_index = try allocator.alloc(u32, constrained_mask.len);
        errdefer allocator.free(reduced_index);
        @memset(reduced_index, std.math.maxInt(u32));

        var free_count: u32 = 0;
        for (owned_mask, 0..) |is_constrained, idx| {
            if (is_constrained) continue;
            reduced_index[idx] = free_count;
            free_count += 1;
        }

        const free_indices = try allocator.alloc(u32, free_count);
        errdefer allocator.free(free_indices);
        for (owned_mask, 0..) |is_constrained, idx| {
            if (is_constrained) continue;
            free_indices[reduced_index[idx]] = @intCast(idx);
        }

        return .{
            .constrained_mask = owned_mask,
            .reduced_index = reduced_index,
            .free_indices = free_indices,
        };
    }

    pub fn initBoundary(
        comptime MeshType: type,
        allocator: std.mem.Allocator,
        mesh: *const MeshType,
        comptime k: comptime_int,
    ) !EliminationMap {
        const mask = try mesh.boundary_mask(allocator, k);
        defer allocator.free(mask);
        return init(allocator, mask);
    }

    pub fn deinit(self: *EliminationMap, allocator: std.mem.Allocator) void {
        allocator.free(self.free_indices);
        allocator.free(self.reduced_index);
        allocator.free(self.constrained_mask);
    }

    pub fn fullCount(self: EliminationMap) usize {
        return self.constrained_mask.len;
    }

    pub fn freeCount(self: EliminationMap) usize {
        return self.free_indices.len;
    }

    pub fn constrainedMask(self: EliminationMap) []const bool {
        return self.constrained_mask;
    }

    pub fn projectFree(self: EliminationMap, full_values: []const f64, reduced_values: []f64) void {
        std.debug.assert(full_values.len == self.fullCount());
        std.debug.assert(reduced_values.len == self.freeCount());
        for (self.free_indices, 0..) |full_idx, reduced_idx| {
            reduced_values[reduced_idx] = full_values[full_idx];
        }
    }

    pub fn liftFree(
        self: EliminationMap,
        reduced_values: []const f64,
        constrained_values: []const f64,
        full_values: []f64,
    ) void {
        std.debug.assert(reduced_values.len == self.freeCount());
        std.debug.assert(constrained_values.len == self.fullCount());
        std.debug.assert(full_values.len == self.fullCount());
        for (self.constrained_mask, 0..) |is_constrained, idx| {
            full_values[idx] = if (is_constrained)
                constrained_values[idx]
            else
                reduced_values[self.reduced_index[idx]];
        }
    }
};

pub const EliminatedLinearSystem = struct {
    elimination_map: EliminationMap,
    boundary_coupling: sparse.CsrMatrix(f64),
    full_rhs: []f64,
    reduced_system: ?LinearSystem,

    pub fn init(
        allocator: std.mem.Allocator,
        full_matrix: sparse.CsrMatrix(f64),
        elimination_map: EliminationMap,
        config: SolveConfig,
    ) !EliminatedLinearSystem {
        std.debug.assert(full_matrix.n_rows == full_matrix.n_cols);
        std.debug.assert(full_matrix.n_rows == elimination_map.fullCount());

        const full_rhs = try allocator.alloc(f64, full_matrix.n_rows);
        errdefer allocator.free(full_rhs);
        @memset(full_rhs, 0.0);

        const reduced_parts = try buildReducedParts(
            allocator,
            full_matrix,
            elimination_map,
            config,
        );
        errdefer {
            reduced_parts.boundary_coupling.deinit(allocator);
            if (reduced_parts.reduced_system) |*system| system.deinit(allocator);
        }

        return .{
            .elimination_map = elimination_map,
            .boundary_coupling = reduced_parts.boundary_coupling,
            .full_rhs = full_rhs,
            .reduced_system = reduced_parts.reduced_system,
        };
    }

    pub fn deinit(self: *EliminatedLinearSystem, allocator: std.mem.Allocator) void {
        if (self.reduced_system) |*system| system.deinit(allocator);
        allocator.free(self.full_rhs);
        self.boundary_coupling.deinit(allocator);
        self.elimination_map.deinit(allocator);
    }

    pub fn fullRhsValues(self: *EliminatedLinearSystem) []f64 {
        return self.full_rhs;
    }

    pub fn seedSolutionFromFull(self: *EliminatedLinearSystem, full_values: []const f64) void {
        std.debug.assert(full_values.len == self.elimination_map.fullCount());
        if (self.reduced_system) |*system| {
            self.elimination_map.projectFree(full_values, system.solutionValues());
        }
    }

    pub fn solveWithConstrainedValues(
        self: *EliminatedLinearSystem,
        constrained_values: []const f64,
        full_solution: []f64,
    ) !cg.SolveResult {
        std.debug.assert(constrained_values.len == self.elimination_map.fullCount());
        std.debug.assert(full_solution.len == self.elimination_map.fullCount());

        if (self.reduced_system) |*system| {
            const reduced_rhs = system.rhsValues();
            for (self.elimination_map.free_indices, 0..) |full_row_idx, reduced_row_idx| {
                var rhs_value = self.full_rhs[full_row_idx];
                const row = self.boundary_coupling.row(@intCast(reduced_row_idx));
                for (row.cols, row.vals) |col_idx, value| {
                    rhs_value -= value * constrained_values[col_idx];
                }
                reduced_rhs[reduced_row_idx] = rhs_value;
            }

            const result = try system.solve();
            self.elimination_map.liftFree(system.solutionValues(), constrained_values, full_solution);
            return result;
        }

        @memcpy(full_solution, constrained_values);
        return .{ .iterations = 0, .relative_residual = 0.0, .converged = true };
    }

    pub fn solveHomogeneous(
        self: *EliminatedLinearSystem,
        full_solution: []f64,
    ) !cg.SolveResult {
        std.debug.assert(full_solution.len == self.elimination_map.fullCount());

        if (self.reduced_system) |*system| {
            self.elimination_map.projectFree(self.full_rhs, system.rhsValues());
            const result = try system.solve();
            for (self.elimination_map.constrained_mask, 0..) |is_constrained, idx| {
                full_solution[idx] = if (is_constrained)
                    0.0
                else
                    system.solutionValues()[self.elimination_map.reduced_index[idx]];
            }
            return result;
        }

        @memset(full_solution, 0.0);
        return .{ .iterations = 0, .relative_residual = 0.0, .converged = true };
    }
};

fn buildReducedParts(
    allocator: std.mem.Allocator,
    full_matrix: sparse.CsrMatrix(f64),
    elimination_map: EliminationMap,
    config: SolveConfig,
) !struct {
    boundary_coupling: sparse.CsrMatrix(f64),
    reduced_system: ?LinearSystem,
} {
    var boundary_triplets = sparse.TripletAssembler(f64).init(
        @intCast(elimination_map.freeCount()),
        full_matrix.n_cols,
    );
    defer boundary_triplets.deinit(allocator);

    if (elimination_map.freeCount() == 0) {
        return .{
            .boundary_coupling = try boundary_triplets.build(allocator),
            .reduced_system = null,
        };
    }

    var reduced_triplets = sparse.TripletAssembler(f64).init(
        @intCast(elimination_map.freeCount()),
        @intCast(elimination_map.freeCount()),
    );
    defer reduced_triplets.deinit(allocator);

    for (0..full_matrix.n_rows) |row_idx_usize| {
        if (elimination_map.constrained_mask[row_idx_usize]) continue;

        const reduced_row = elimination_map.reduced_index[row_idx_usize];
        const row = full_matrix.row(@intCast(row_idx_usize));
        for (row.cols, row.vals) |col_idx, value| {
            if (elimination_map.constrained_mask[col_idx]) {
                try boundary_triplets.addEntry(allocator, reduced_row, col_idx, value);
            } else {
                try reduced_triplets.addEntry(allocator, reduced_row, elimination_map.reduced_index[col_idx], value);
            }
        }
    }

    var reduced_matrix = try reduced_triplets.build(allocator);
    errdefer reduced_matrix.deinit(allocator);

    var boundary_coupling = try boundary_triplets.build(allocator);
    errdefer boundary_coupling.deinit(allocator);

    return .{
        .boundary_coupling = boundary_coupling,
        .reduced_system = try LinearSystem.init(allocator, reduced_matrix, config),
    };
}

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

test "linear system computes matrix diagonal and reuses owned buffers across solves" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 2, 2, 2);
    matrix.row_ptr[0] = 0;
    matrix.row_ptr[1] = 1;
    matrix.row_ptr[2] = 2;
    matrix.col_idx[0] = 0;
    matrix.col_idx[1] = 1;
    matrix.values[0] = 2.0;
    matrix.values[1] = 4.0;

    var system = try LinearSystem.init(allocator, matrix, .{});
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

test "eliminated linear system scatters constrained values and corrects reduced rhs" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 3, 3, 7);
    defer matrix.deinit(allocator);
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

    const elimination_map = try EliminationMap.init(allocator, &[_]bool{ true, false, true });
    var system = try EliminatedLinearSystem.init(allocator, matrix, elimination_map, .{});
    defer system.deinit(allocator);

    system.fullRhsValues()[0] = 0.0;
    system.fullRhsValues()[1] = 0.0;
    system.fullRhsValues()[2] = 0.0;

    const constrained_values = [_]f64{ 10.0, 0.0, 20.0 };
    var solution = [_]f64{ 0.0, 0.0, 0.0 };

    const result = try system.solveWithConstrainedValues(&constrained_values, &solution);
    try testing.expect(result.converged);
    try testing.expectApproxEqAbs(@as(f64, 10.0), solution[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 15.0), solution[1], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 20.0), solution[2], 1e-12);
}

test "eliminated linear system reuses full rhs buffer and seeds from full state" {
    const allocator = testing.allocator;

    var matrix = try sparse.CsrMatrix(f64).init(allocator, 3, 3, 3);
    defer matrix.deinit(allocator);
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

    const elimination_map = try EliminationMap.init(allocator, &[_]bool{ true, false, false });
    var system = try EliminatedLinearSystem.init(allocator, matrix, elimination_map, .{});
    defer system.deinit(allocator);

    const rhs_ptr = system.fullRhsValues().ptr;
    system.seedSolutionFromFull(&[_]f64{ 100.0, 10.0, 18.0 });
    system.fullRhsValues()[0] = 0.0;
    system.fullRhsValues()[1] = 10.0;
    system.fullRhsValues()[2] = 18.0;

    var solution = [_]f64{ 0.0, 0.0, 0.0 };
    const result = try system.solveHomogeneous(&solution);
    try testing.expect(result.converged);
    try testing.expect(rhs_ptr == system.fullRhsValues().ptr);
    try testing.expectApproxEqAbs(@as(f64, 0.0), solution[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 2.0), solution[1], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 3.0), solution[2], 1e-12);
}

test {
    testing.refAllDeclsRecursive(@This());
}
