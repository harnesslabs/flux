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
        _ = allocator;
        _ = matrix;
        _ = config;
        @panic("not yet implemented");
    }

    pub fn deinit(self: *AssembledImplicitSystem, allocator: std.mem.Allocator) void {
        _ = self;
        _ = allocator;
        @panic("not yet implemented");
    }

    pub fn rhsValues(self: *AssembledImplicitSystem) []f64 {
        _ = self;
        @panic("not yet implemented");
    }

    pub fn solutionValues(self: *AssembledImplicitSystem) []f64 {
        _ = self;
        @panic("not yet implemented");
    }

    pub fn seedSolution(self: *AssembledImplicitSystem, values: []const f64) void {
        _ = self;
        _ = values;
        @panic("not yet implemented");
    }

    pub fn solve(self: *AssembledImplicitSystem) !cg.SolveResult {
        _ = self;
        @panic("not yet implemented");
    }
};

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

test {
    testing.refAllDeclsRecursive(@This());
}
