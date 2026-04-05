//! Poisson solve helpers for primal 0-forms with Dirichlet boundary conditions.
//!
//! The discrete scalar Laplacian satisfies `Δ₀ = ★₀⁻¹ S`, where
//! `S = D₀ᵀ M₁ D₀` is symmetric positive semidefinite. Solving `Δ₀u = f`
//! is therefore equivalent to solving `S u = ★₀ f`, then enforcing Dirichlet
//! values by symmetric row/column elimination so the reduced system remains SPD.

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");
const obj = @import("../io/obj.zig");
const operator_context_mod = @import("context.zig");
const conjugate_gradient = @import("../math/cg.zig");

comptime {
    @setEvalBranchQuota(20000);
}

pub const SolveConfig = struct {
    tolerance_relative: f64 = 1e-10,
    iteration_limit: u32 = 4000,
};

pub const SolveResult = struct {
    solution: []f64,
    cg_result: conjugate_gradient.SolveResult,

    pub fn deinit(self: *SolveResult, allocator: std.mem.Allocator) void {
        allocator.free(self.solution);
    }
};

pub fn solve_zero_form_dirichlet(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    operator_context: anytype,
    forcing_values: []const f64,
    boundary_values: []const f64,
    config: SolveConfig,
) !SolveResult {
    _ = MeshType;
    _ = allocator;
    _ = operator_context;
    _ = forcing_values;
    _ = boundary_values;
    _ = config;
    return error.NotYetImplemented;
}

fn exact_solution_2d(coords: [2]f64) f64 {
    return std.math.sin(std.math.pi * coords[0]) * std.math.sin(std.math.pi * coords[1]);
}

fn rhs_2d(coords: [2]f64) f64 {
    return 2.0 * std.math.pi * std.math.pi * exact_solution_2d(coords);
}

fn exact_solution_3d(coords: [3]f64) f64 {
    return std.math.sin(std.math.pi * coords[0]) *
        std.math.sin(std.math.pi * coords[1]) *
        std.math.sin(std.math.pi * coords[2]);
}

fn rhs_3d(coords: [3]f64) f64 {
    return 3.0 * std.math.pi * std.math.pi * exact_solution_3d(coords);
}

fn weighted_l2_error(mesh: anytype, approx: []const f64, exact: []const f64) f64 {
    std.debug.assert(approx.len == exact.len);
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    var error_sq: f64 = 0.0;
    var measure: f64 = 0.0;
    for (approx, exact, dual_volumes) |uh, u, dual_volume| {
        const diff = uh - u;
        error_sq += diff * diff * dual_volume;
        measure += dual_volume;
    }
    return @sqrt(error_sq / measure);
}

fn expect_second_order(errors: []const f64) !void {
    std.debug.assert(errors.len >= 2);
    for (0..errors.len - 1) |idx| {
        const rate = std.math.log(f64, 2.0, errors[idx] / errors[idx + 1]);
        try testing.expect(rate > 1.75);
    }
}

fn fill_vertex_samples_2d(mesh: *const topology.Mesh(2, 2), forcing: []f64, boundary: []f64, exact: []f64) void {
    const coords = mesh.vertices.slice().items(.coords);
    for (coords, forcing, boundary, exact) |coord, *rhs_value, *boundary_value, *exact_value| {
        const xy = [2]f64{ coord[0], coord[1] };
        exact_value.* = exact_solution_2d(xy);
        rhs_value.* = rhs_2d(xy);
        boundary_value.* = exact_value.*;
    }
}

fn fill_vertex_samples_surface(mesh: *const topology.Mesh(3, 2), forcing: []f64, boundary: []f64, exact: []f64) void {
    const coords = mesh.vertices.slice().items(.coords);
    for (coords, forcing, boundary, exact) |coord, *rhs_value, *boundary_value, *exact_value| {
        const xy = [2]f64{ coord[0], coord[1] };
        exact_value.* = exact_solution_2d(xy);
        rhs_value.* = rhs_2d(xy);
        boundary_value.* = exact_value.*;
    }
}

fn fill_vertex_samples_3d(mesh: *const topology.Mesh(3, 3), forcing: []f64, boundary: []f64, exact: []f64) void {
    const coords = mesh.vertices.slice().items(.coords);
    for (coords, forcing, boundary, exact) |coord, *rhs_value, *boundary_value, *exact_value| {
        exact_value.* = exact_solution_3d(coord);
        rhs_value.* = rhs_3d(coord);
        boundary_value.* = exact_value.*;
    }
}

fn make_planar_quad_obj(allocator: std.mem.Allocator, nx: u32, ny: u32) ![]u8 {
    var output = std.ArrayListUnmanaged(u8){};
    errdefer output.deinit(allocator);

    for (0..nx + 1) |i_u| {
        const x = @as(f64, @floatFromInt(i_u)) / @as(f64, @floatFromInt(nx));
        for (0..ny + 1) |j_u| {
            const y = @as(f64, @floatFromInt(j_u)) / @as(f64, @floatFromInt(ny));
            try output.writer(allocator).print("v {d} {d} 0\n", .{ x, y });
        }
    }

    for (0..nx) |i_u| {
        for (0..ny) |j_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            const row_stride = ny + 1;
            const sw = i * row_stride + j + 1;
            const se = (i + 1) * row_stride + j + 1;
            const ne = (i + 1) * row_stride + (j + 1) + 1;
            const nw = i * row_stride + (j + 1) + 1;
            try output.writer(allocator).print("f {d} {d} {d} {d}\n", .{ sw, se, ne, nw });
        }
    }

    return output.toOwnedSlice(allocator);
}

test "Poisson solve converges at second order on 2D uniform grids" {
    const allocator = testing.allocator;
    const Mesh2D = topology.Mesh(2, 2);
    const levels = [_]u32{ 4, 8, 16 };
    var errors: [levels.len]f64 = undefined;

    for (levels, 0..) |n, level_idx| {
        var mesh = try Mesh2D.uniform_grid(allocator, n, n, 1.0, 1.0);
        defer mesh.deinit(allocator);

        var operator_context = try operator_context_mod.OperatorContext(Mesh2D).init(allocator, &mesh);
        defer operator_context.deinit();

        const vertex_count = mesh.num_vertices();
        const forcing = try allocator.alloc(f64, vertex_count);
        defer allocator.free(forcing);
        const boundary = try allocator.alloc(f64, vertex_count);
        defer allocator.free(boundary);
        const exact = try allocator.alloc(f64, vertex_count);
        defer allocator.free(exact);
        fill_vertex_samples_2d(&mesh, forcing, boundary, exact);

        var solve = try solve_zero_form_dirichlet(Mesh2D, allocator, operator_context, forcing, boundary, .{});
        defer solve.deinit(allocator);

        errors[level_idx] = weighted_l2_error(&mesh, solve.solution, exact);
    }

    try expect_second_order(&errors);
}

test "Poisson solve converges at second order on tetrahedral grids" {
    const allocator = testing.allocator;
    const Mesh3D = topology.Mesh(3, 3);
    const levels = [_]u32{ 2, 4, 8 };
    var errors: [levels.len]f64 = undefined;

    for (levels, 0..) |n, level_idx| {
        var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, n, n, n, 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        var operator_context = try operator_context_mod.OperatorContext(Mesh3D).init(allocator, &mesh);
        defer operator_context.deinit();

        const vertex_count = mesh.num_vertices();
        const forcing = try allocator.alloc(f64, vertex_count);
        defer allocator.free(forcing);
        const boundary = try allocator.alloc(f64, vertex_count);
        defer allocator.free(boundary);
        const exact = try allocator.alloc(f64, vertex_count);
        defer allocator.free(exact);
        fill_vertex_samples_3d(&mesh, forcing, boundary, exact);

        var solve = try solve_zero_form_dirichlet(Mesh3D, allocator, operator_context, forcing, boundary, .{});
        defer solve.deinit(allocator);

        errors[level_idx] = weighted_l2_error(&mesh, solve.solution, exact);
    }

    try expect_second_order(&errors);
}

test "Poisson solve converges on imported quad OBJ surface embedded in R3" {
    const allocator = testing.allocator;
    const MeshSurface = topology.Mesh(3, 2);
    const levels = [_]u32{ 4, 8, 16 };
    var errors: [levels.len]f64 = undefined;

    for (levels, 0..) |n, level_idx| {
        const obj_bytes = try make_planar_quad_obj(allocator, n, n);
        defer allocator.free(obj_bytes);

        var raw_mesh = try obj.parse_obj(allocator, obj_bytes);
        defer raw_mesh.deinit(allocator);
        var triangulated = try obj.triangulate(allocator, &raw_mesh);
        defer triangulated.deinit(allocator);

        var mesh = try MeshSurface.from_triangles(allocator, triangulated.vertex_coords, triangulated.triangles);
        defer mesh.deinit(allocator);

        var operator_context = try operator_context_mod.OperatorContext(MeshSurface).init(allocator, &mesh);
        defer operator_context.deinit();

        const vertex_count = mesh.num_vertices();
        const forcing = try allocator.alloc(f64, vertex_count);
        defer allocator.free(forcing);
        const boundary = try allocator.alloc(f64, vertex_count);
        defer allocator.free(boundary);
        const exact = try allocator.alloc(f64, vertex_count);
        defer allocator.free(exact);
        fill_vertex_samples_surface(&mesh, forcing, boundary, exact);

        var solve = try solve_zero_form_dirichlet(MeshSurface, allocator, operator_context, forcing, boundary, .{});
        defer solve.deinit(allocator);

        errors[level_idx] = weighted_l2_error(&mesh, solve.solution, exact);
    }

    try expect_second_order(&errors);
}
