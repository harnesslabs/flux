//! DEC strong-operator surfaces for the Laplace-de Rham operator.
//!
//! This module owns the strong-action path `Δ = dδ + δd`. When a supported
//! FEEC weak stiffness matrix exists, the assembled DEC operator caches it as a
//! solve-oriented auxiliary object, but the primary semantic here remains
//! strong operator application on primal cochains.

const std = @import("std");
const cochain = @import("../../forms/cochain.zig");
const exterior_derivative = @import("../exterior_derivative.zig");
const hodge_star = @import("../hodge_star.zig");
const sparse = @import("../../math/sparse.zig");
const feec = @import("feec.zig");

fn validatePrimalCochainInput(comptime InputType: type) void {
    if (!@hasDecl(InputType, "duality")) {
        @compileError("laplacian requires a Cochain type");
    }
    if (InputType.duality != cochain.Primal) {
        @compileError("laplacian expects a primal cochain");
    }
}

/// Stored DEC Laplacian operator specialized to a primal k-cochain type.
///
/// The operator applies the strong action `Δₖ`. For supported FEEC weak-form
/// cases it also keeps the corresponding stiffness matrix `Sₖ` as auxiliary
/// assembled state for solver and consistency-check paths.
pub fn AssembledLaplacian(comptime InputType: type) type {
    comptime validatePrimalCochainInput(InputType);

    return struct {
        const Self = @This();
        const MeshType = InputType.MeshT;
        const k = InputType.degree;

        mesh: *const MeshType,
        stiffness: sparse.CsrMatrix(f64),
        left_scaling: []f64,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.left_scaling);
            self.stiffness.deinit(allocator);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            std.debug.assert(input.mesh == self.mesh);

            var output = try InputType.init(allocator, input.mesh);
            errdefer output.deinit(allocator);

            switch (k) {
                0 => {
                    sparse.spmv(self.stiffness, input.values, output.values);
                    for (output.values, self.left_scaling) |*out, scale| {
                        out.* *= scale;
                    }
                    return output;
                },
                else => {
                    output.deinit(allocator);
                    return laplacian_composed(allocator, input);
                },
            }
        }
    };
}

/// Assemble a stored DEC Laplacian operator for a fixed mesh and primal degree.
pub fn assemble_for_degree(
    comptime MeshType: type,
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledLaplacian(cochain.Cochain(MeshType, k, cochain.Primal)) {
    const stiffness_supported = comptime feec.supportsExplicitStiffness(MeshType, k);
    var stiffness = switch (stiffness_supported) {
        true => try feec.assemble_stiffness(k, allocator, mesh),
        false => try sparse.CsrMatrix(f64).init(allocator, 0, 0, 0),
    };
    errdefer stiffness.deinit(allocator);

    const left_scaling = switch (k) {
        0 => try assemble_zero_form_star_inverse_diag(allocator, mesh),
        else => try allocator.alloc(f64, 0),
    };
    errdefer allocator.free(left_scaling);

    return .{
        .mesh = mesh,
        .stiffness = stiffness,
        .left_scaling = left_scaling,
    };
}

fn assemble_zero_form_star_inverse_diag(
    allocator: std.mem.Allocator,
    mesh: anytype,
) ![]f64 {
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    const diagonal = try allocator.alloc(f64, dual_volumes.len);
    errdefer allocator.free(diagonal);

    for (diagonal, dual_volumes) |*out, dual_volume| {
        std.debug.assert(dual_volume != 0.0);
        out.* = 1.0 / dual_volume;
    }

    return diagonal;
}

pub fn laplacian_composed(
    allocator: std.mem.Allocator,
    input: anytype,
) !@TypeOf(input) {
    const InputType = @TypeOf(input);
    comptime validatePrimalCochainInput(InputType);

    const MeshType = InputType.MeshT;
    const k = InputType.degree;
    const topological_dimension = MeshType.topological_dimension;

    var result = try InputType.init(allocator, input.mesh);
    errdefer result.deinit(allocator);

    const bk1_cols = if (k < topological_dimension) input.mesh.boundary(k + 1).n_cols else 0;
    const bk_cols = if (k > 0) input.mesh.boundary(k).n_cols else 0;
    const workspace_len = bk1_cols + 2 * bk_cols;

    const workspace = try allocator.alloc(f64, workspace_len);
    defer allocator.free(workspace);

    if (k < topological_dimension) {
        var d_omega = try exterior_derivative.exterior_derivative(allocator, input);
        defer d_omega.deinit(allocator);

        var star_d = try hodge_star.hodge_star(allocator, d_omega);
        defer star_d.deinit(allocator);

        const bk1 = input.mesh.boundary(k + 1);
        const temp = workspace[0..bk1_cols];
        @memset(temp, 0);
        bk1.transpose_multiply(star_d.values, temp);

        try hodge_star.apply_inverse_raw(allocator, MeshType, k, input.mesh, temp, result.values);
    }

    if (k > 0) {
        var star_omega = try hodge_star.hodge_star(allocator, input);
        defer star_omega.deinit(allocator);

        const bk = input.mesh.boundary(k);
        const temp_km1 = workspace[bk1_cols .. bk1_cols + bk_cols];
        @memset(temp_km1, 0);
        bk.transpose_multiply(star_omega.values, temp_km1);

        const codiff_vals = workspace[bk1_cols + bk_cols ..];
        try hodge_star.apply_inverse_raw(allocator, MeshType, k - 1, input.mesh, temp_km1, codiff_vals);

        for (0..bk.n_rows) |row_idx| {
            const r = bk.row(@intCast(row_idx));
            var sum: f64 = 0;
            for (r.cols, 0..) |col, entry_idx| {
                const sign = r.sign(entry_idx);
                sum += @as(f64, @floatFromInt(sign)) * codiff_vals[col];
            }
            result.values[row_idx] += sum;
        }
    }

    return result;
}
