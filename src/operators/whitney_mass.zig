//! Whitney mass matrices for Galerkin Hodge stars on interior degrees.
//!
//! On barycentric duals, the diagonal DEC ratio `dual_measure / primal_measure`
//! is a good preconditioner but not the correct Galerkin Hodge star away from
//! orthogonal primal/dual pairs. The Whitney mass matrix replaces that diagonal
//! approximation with the FEEC inner product
//!
//!   M_k(σᵢ, σⱼ) = ∫_Ω W_{σᵢ} · W_{σⱼ} dV
//!
//! for interior degrees `0 < k < n`.
//!
//! References:
//!   - Bossavit, "Whitney forms: a class of finite elements for three-dimensional
//!     computations in electromagnetism" (1988)
//!   - Arnold, Falk, Winther, "Finite Element Exterior Calculus" (2006), §5
//!   - Hirani, "Discrete Exterior Calculus" (2003), §2.4

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");
const sparse = @import("../math/sparse.zig");

pub fn assemble_whitney_mass(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const MeshType = @TypeOf(mesh.*);
    const n = MeshType.topological_dimension;

    comptime {
        if (k <= 0 or k >= n) {
            @compileError("Whitney mass is only defined for interior degrees 0 < k < n");
        }
        if (n > 3) {
            @compileError("Whitney mass assembly is currently implemented for topological_dimension <= 3");
        }
    }

    const simplex_count = mesh.num_cells(k);
    const top_simplex_count = mesh.num_cells(n);
    const simplex_vertices = mesh.simplices(k).items(.vertices);
    const top_simplex_vertices = mesh.simplices(n).items(.vertices);
    const top_simplex_volumes = mesh.simplices(n).items(.volume);
    const coords = mesh.vertices.slice().items(.coords);
    const local_faces = localFaces(n, k);

    var simplex_index = std.AutoHashMap([k + 1]u32, u32).init(allocator);
    defer simplex_index.deinit();
    for (simplex_vertices, 0..) |vertices, global_idx| {
        try simplex_index.put(vertices, @intCast(global_idx));
    }

    var assembler = sparse.TripletAssembler(f64).init(simplex_count, simplex_count);
    defer assembler.deinit(allocator);

    for (0..top_simplex_count) |top_idx| {
        const top_vertices = top_simplex_vertices[top_idx];

        var top_coords: [n + 1][MeshType.embedding_dimension]f64 = undefined;
        for (0..n + 1) |local_vertex_idx| {
            top_coords[local_vertex_idx] = coords[top_vertices[local_vertex_idx]];
        }

        const gradients = barycentricGradients(MeshType.embedding_dimension, n, top_coords);
        const local_mass = localWhitneyMass(MeshType.embedding_dimension, n, k, gradients, top_simplex_volumes[top_idx]);

        var global_indices: [local_faces.len]u32 = undefined;
        var orientation_signs: [local_faces.len]i8 = undefined;
        for (local_faces, 0..) |local_face, local_face_idx| {
            const oriented_vertices = liftLocalFaceVertices(k, n, top_vertices, local_face);
            const canonical_key = canonicalizeVertices(k + 1, oriented_vertices);
            const global_idx = simplex_index.get(canonical_key).?;
            global_indices[local_face_idx] = global_idx;
            orientation_signs[local_face_idx] = orientationSign(k + 1, oriented_vertices, simplex_vertices[global_idx]);
        }

        for (0..local_faces.len) |i| {
            const sign_i: f64 = @floatFromInt(orientation_signs[i]);
            for (0..local_faces.len) |j| {
                const sign_j: f64 = @floatFromInt(orientation_signs[j]);
                try assembler.addEntry(allocator, global_indices[i], global_indices[j], sign_i * sign_j * local_mass[i][j]);
            }
        }
    }

    return assembler.build(allocator);
}

pub fn assemble_whitney_mass_1(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    return assemble_whitney_mass(1, allocator, mesh);
}

pub fn assemble_whitney_preconditioner(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
) ![]f64 {
    const MeshType = @TypeOf(mesh.*);
    const n = MeshType.topological_dimension;

    comptime {
        if (k <= 0 or k >= n) {
            @compileError("Whitney preconditioner is only defined for interior degrees 0 < k < n");
        }
    }

    const preconditioner = try allocator.alloc(f64, mesh.num_cells(k));
    errdefer allocator.free(preconditioner);

    if (k == 1) {
        const primal_measures = mesh.simplices(1).items(.volume);
        for (preconditioner, primal_measures, mesh.dual_edge_volumes) |*out, primal_measure, dual_measure| {
            std.debug.assert(primal_measure != 0.0);
            out.* = dual_measure / primal_measure;
        }
        return preconditioner;
    }

    if (n == 3 and k == 2) {
        const dual_face_lengths = try allocator.alloc(f64, mesh.num_faces());
        defer allocator.free(dual_face_lengths);
        accumulateDualFaceLengths(mesh, dual_face_lengths);

        const primal_measures = mesh.simplices(2).items(.volume);
        for (preconditioner, primal_measures, dual_face_lengths) |*out, primal_measure, dual_measure| {
            std.debug.assert(primal_measure != 0.0);
            out.* = dual_measure / primal_measure;
        }
        return preconditioner;
    }

    unreachable;
}

fn vecDot(comptime dimension: usize, a: [dimension]f64, b: [dimension]f64) f64 {
    var sum: f64 = 0;
    for (0..dimension) |idx| sum += a[idx] * b[idx];
    return sum;
}

fn vecSub(comptime dimension: usize, a: [dimension]f64, b: [dimension]f64) [dimension]f64 {
    var result: [dimension]f64 = undefined;
    for (0..dimension) |idx| result[idx] = a[idx] - b[idx];
    return result;
}

fn factorial(comptime n: comptime_int) comptime_int {
    var result: comptime_int = 1;
    inline for (1..n + 1) |i| result *= i;
    return result;
}

fn choose(comptime n: comptime_int, comptime k: comptime_int) comptime_int {
    if (k < 0 or k > n) return 0;
    if (k == 0 or k == n) return 1;
    var numerator: comptime_int = 1;
    var denominator: comptime_int = 1;
    const k_small = if (k < n - k) k else n - k;
    inline for (0..k_small) |i| {
        numerator *= (n - i);
        denominator *= (i + 1);
    }
    return @divExact(numerator, denominator);
}

fn localFaces(comptime n: comptime_int, comptime k: comptime_int) [choose(n + 1, k + 1)][k + 1]u8 {
    return switch (n) {
        2 => switch (k) {
            1 => .{
                .{ 0, 1 },
                .{ 0, 2 },
                .{ 1, 2 },
            },
            else => @compileError("unsupported local face degree for 2-simplex"),
        },
        3 => switch (k) {
            1 => .{
                .{ 0, 1 },
                .{ 0, 2 },
                .{ 0, 3 },
                .{ 1, 2 },
                .{ 1, 3 },
                .{ 2, 3 },
            },
            2 => .{
                .{ 0, 1, 2 },
                .{ 0, 1, 3 },
                .{ 0, 2, 3 },
                .{ 1, 2, 3 },
            },
            else => @compileError("unsupported local face degree for 3-simplex"),
        },
        else => @compileError("local face enumeration is only implemented for topological_dimension <= 3"),
    };
}

fn liftLocalFaceVertices(
    comptime k: comptime_int,
    comptime n: comptime_int,
    top_vertices: [n + 1]u32,
    local_face: [k + 1]u8,
) [k + 1]u32 {
    var result: [k + 1]u32 = undefined;
    for (local_face, 0..) |local_vertex_idx, write_idx| {
        result[write_idx] = top_vertices[local_vertex_idx];
    }
    return result;
}

fn canonicalizeVertices(comptime len: comptime_int, vertices: [len]u32) [len]u32 {
    var result = vertices;
    inline for (1..len) |i| {
        var j = i;
        while (j > 0 and result[j - 1] > result[j]) : (j -= 1) {
            std.mem.swap(u32, &result[j - 1], &result[j]);
        }
    }
    return result;
}

fn orientationSign(comptime len: comptime_int, oriented_vertices: [len]u32, global_vertices: [len]u32) i8 {
    var permutation: [len]u8 = undefined;
    var used = [_]bool{false} ** len;
    for (oriented_vertices, 0..) |vertex, oriented_idx| {
        var found = false;
        for (global_vertices, 0..) |global_vertex, global_idx| {
            if (used[global_idx]) continue;
            if (vertex != global_vertex) continue;
            permutation[oriented_idx] = @intCast(global_idx);
            used[global_idx] = true;
            found = true;
            break;
        }
        std.debug.assert(found);
    }

    var inversion_count: u32 = 0;
    inline for (0..len) |i| {
        inline for (i + 1..len) |j| {
            if (permutation[i] > permutation[j]) inversion_count += 1;
        }
    }
    return if (inversion_count % 2 == 0) 1 else -1;
}

fn barycentricGradients(
    comptime embedding_dimension: usize,
    comptime n: comptime_int,
    simplex_vertices: [n + 1][embedding_dimension]f64,
) [n + 1][embedding_dimension]f64 {
    var edge_vectors: [n][embedding_dimension]f64 = undefined;
    inline for (1..n + 1) |local_vertex_idx| {
        edge_vectors[local_vertex_idx - 1] = vecSub(embedding_dimension, simplex_vertices[local_vertex_idx], simplex_vertices[0]);
    }

    var gram: [n][n]f64 = undefined;
    inline for (0..n) |i| {
        inline for (0..n) |j| {
            gram[i][j] = vecDot(embedding_dimension, edge_vectors[i], edge_vectors[j]);
        }
    }
    const gram_inverse = invertSmallMatrix(n, gram);

    var gradients: [n + 1][embedding_dimension]f64 = undefined;
    gradients[0] = @splat(0.0);
    inline for (1..n + 1) |local_vertex_idx| {
        var gradient: [embedding_dimension]f64 = @splat(0.0);
        inline for (0..n) |edge_idx| {
            const coefficient = gram_inverse[edge_idx][local_vertex_idx - 1];
            inline for (0..embedding_dimension) |dim_idx| {
                gradient[dim_idx] += edge_vectors[edge_idx][dim_idx] * coefficient;
            }
        }
        gradients[local_vertex_idx] = gradient;
        inline for (0..embedding_dimension) |dim_idx| {
            gradients[0][dim_idx] -= gradient[dim_idx];
        }
    }

    return gradients;
}

fn localWhitneyMass(
    comptime embedding_dimension: usize,
    comptime n: comptime_int,
    comptime k: comptime_int,
    gradients: [n + 1][embedding_dimension]f64,
    simplex_volume: f64,
) [choose(n + 1, k + 1)][choose(n + 1, k + 1)]f64 {
    const local_faces = localFaces(n, k);
    var local_mass: [local_faces.len][local_faces.len]f64 = undefined;
    const lambda_integral_scale = simplex_volume / @as(f64, @floatFromInt((n + 1) * (n + 2)));
    const whitney_scale = @as(f64, @floatFromInt(factorial(k) * factorial(k)));

    for (local_faces, 0..) |left_face, left_idx| {
        for (local_faces, 0..) |right_face, right_idx| {
            var entry: f64 = 0.0;
            inline for (0..k + 1) |left_omit| {
                const left_basis = omitIndex(k + 1, left_face, left_omit);
                const left_lambda = left_face[left_omit];
                inline for (0..k + 1) |right_omit| {
                    const right_basis = omitIndex(k + 1, right_face, right_omit);
                    const right_lambda = right_face[right_omit];
                    const sign: f64 = if ((left_omit + right_omit) % 2 == 0) 1.0 else -1.0;
                    const lambda_inner: f64 = if (left_lambda == right_lambda) 2.0 else 1.0;
                    entry += sign * lambda_inner *
                        wedgeInnerProduct(embedding_dimension, k, gradients[0..], left_basis, right_basis);
                }
            }
            local_mass[left_idx][right_idx] = whitney_scale * lambda_integral_scale * entry;
        }
    }

    return local_mass;
}

fn omitIndex(comptime len: comptime_int, indices: [len]u8, omit: usize) [len - 1]u8 {
    var result: [len - 1]u8 = undefined;
    var write_idx: usize = 0;
    for (indices, 0..) |index, read_idx| {
        if (read_idx == omit) continue;
        result[write_idx] = index;
        write_idx += 1;
    }
    return result;
}

fn wedgeInnerProduct(
    comptime embedding_dimension: usize,
    comptime k: comptime_int,
    gradients: []const [embedding_dimension]f64,
    left_indices: [k]u8,
    right_indices: [k]u8,
) f64 {
    if (k == 1) {
        return vecDot(embedding_dimension, gradients[left_indices[0]], gradients[right_indices[0]]);
    }

    if (k == 2) {
        const a = vecDot(embedding_dimension, gradients[left_indices[0]], gradients[right_indices[0]]);
        const b = vecDot(embedding_dimension, gradients[left_indices[0]], gradients[right_indices[1]]);
        const c = vecDot(embedding_dimension, gradients[left_indices[1]], gradients[right_indices[0]]);
        const d = vecDot(embedding_dimension, gradients[left_indices[1]], gradients[right_indices[1]]);
        return a * d - b * c;
    }

    @compileError("Whitney wedge inner product is only implemented for k <= 2");
}

fn invertSmallMatrix(comptime n: comptime_int, matrix: [n][n]f64) [n][n]f64 {
    if (n == 2) {
        const det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        std.debug.assert(det != 0.0);
        return .{
            .{ matrix[1][1] / det, -matrix[0][1] / det },
            .{ -matrix[1][0] / det, matrix[0][0] / det },
        };
    }

    if (n == 3) {
        const a = matrix[0][0];
        const b = matrix[0][1];
        const c = matrix[0][2];
        const d = matrix[1][0];
        const e = matrix[1][1];
        const f = matrix[1][2];
        const g = matrix[2][0];
        const h = matrix[2][1];
        const i = matrix[2][2];

        const cofactor00 = e * i - f * h;
        const cofactor01 = -(d * i - f * g);
        const cofactor02 = d * h - e * g;
        const cofactor10 = -(b * i - c * h);
        const cofactor11 = a * i - c * g;
        const cofactor12 = -(a * h - b * g);
        const cofactor20 = b * f - c * e;
        const cofactor21 = -(a * f - c * d);
        const cofactor22 = a * e - b * d;

        const det = a * cofactor00 + b * cofactor01 + c * cofactor02;
        std.debug.assert(det != 0.0);

        return .{
            .{ cofactor00 / det, cofactor10 / det, cofactor20 / det },
            .{ cofactor01 / det, cofactor11 / det, cofactor21 / det },
            .{ cofactor02 / det, cofactor12 / det, cofactor22 / det },
        };
    }

    @compileError("small matrix inversion is only implemented for n = 2 or n = 3");
}

fn accumulateDualFaceLengths(
    mesh: *const topology.Mesh(3, 3),
    dual_face_lengths: []f64,
) void {
    std.debug.assert(dual_face_lengths.len == mesh.num_faces());
    @memset(dual_face_lengths, 0.0);

    const face_barycenters = mesh.simplices(2).items(.barycenter);
    const tet_barycenters = mesh.simplices(3).items(.barycenter);
    for (0..mesh.num_tets()) |tet_idx| {
        const tet_barycenter = tet_barycenters[tet_idx];
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx));
        for (tet_row.cols) |face_idx| {
            dual_face_lengths[face_idx] += euclideanDistance(3, face_barycenters[face_idx], tet_barycenter);
        }
    }
}

fn euclideanDistance(comptime dimension: usize, a: [dimension]f64, b: [dimension]f64) f64 {
    var sum_squared: f64 = 0.0;
    for (0..dimension) |idx| {
        const delta = a[idx] - b[idx];
        sum_squared += delta * delta;
    }
    return @sqrt(sum_squared);
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
const Mesh3D = topology.Mesh(3, 3);

test "Whitney mass matrix is square with n_edges dimension" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try assemble_whitney_mass_1(allocator, &mesh);
    defer mass.deinit(allocator);

    try testing.expectEqual(mesh.num_edges(), mass.n_rows);
    try testing.expectEqual(mesh.num_edges(), mass.n_cols);
}

test "Whitney mass matrix is symmetric" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var mass = try assemble_whitney_mass_1(allocator, &mesh);
    defer mass.deinit(allocator);

    // Check M(i,j) = M(j,i) by comparing Ax with Aᵀx for random x.
    var rng = std.Random.DefaultPrng.init(0xDEC_AA55_00);

    const x = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(y);
    const ax = try allocator.alloc(f64, mass.n_rows);
    defer allocator.free(ax);
    const atx = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(atx);

    for (0..100) |_| {
        for (x) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;
        for (y) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        // Compute ⟨Ax, y⟩ and ⟨x, Ay⟩ — they must be equal for symmetric A.
        sparse.spmv(mass, x, ax);
        @memset(atx, 0.0);

        // ⟨Ax, y⟩
        var dot_axy: f64 = 0;
        for (ax, y) |a, b| dot_axy += a * b;

        // ⟨x, Ay⟩
        sparse.spmv(mass, y, atx);
        var dot_xay: f64 = 0;
        for (x, atx) |a, b| dot_xay += a * b;

        try testing.expectApproxEqRel(dot_axy, dot_xay, 1e-12);
    }
}

test "Whitney mass matrix is positive definite (xᵀMx > 0 for random x)" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try assemble_whitney_mass_1(allocator, &mesh);
    defer mass.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xDEC_AA55_01);

    const x = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(x);
    const mx = try allocator.alloc(f64, mass.n_rows);
    defer allocator.free(mx);

    for (0..1000) |_| {
        for (x) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        sparse.spmv(mass, x, mx);

        var xtmx: f64 = 0;
        for (x, mx) |xi, mxi| xtmx += xi * mxi;

        try testing.expect(xtmx > 0.0);
    }
}

test "Whitney mass matrix has nonzero diagonal and at least one entry per row" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try assemble_whitney_mass_1(allocator, &mesh);
    defer mass.deinit(allocator);

    // Every edge appears in at least one face, so every row has entries.
    // The diagonal (self-inner-product) must be strictly positive.
    for (0..mass.n_rows) |row_idx| {
        const r = mass.row(@intCast(row_idx));
        try testing.expect(r.cols.len > 0);

        // Find diagonal entry.
        var found_diag = false;
        for (r.cols, r.vals) |col, val| {
            if (col == row_idx) {
                try testing.expect(val > 0.0);
                found_diag = true;
            }
        }
        try testing.expect(found_diag);
    }
}

test "Whitney 2-form mass matrix is symmetric positive definite on tetrahedral mesh" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 2, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var mass = try assemble_whitney_mass(2, allocator, &mesh);
    defer mass.deinit(allocator);

    try testing.expectEqual(mesh.num_faces(), mass.n_rows);
    try testing.expectEqual(mesh.num_faces(), mass.n_cols);

    var rng = std.Random.DefaultPrng.init(0xDEC_3D_02);
    const x = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, mass.n_cols);
    defer allocator.free(y);
    const ax = try allocator.alloc(f64, mass.n_rows);
    defer allocator.free(ax);
    const ay = try allocator.alloc(f64, mass.n_rows);
    defer allocator.free(ay);

    for (0..100) |_| {
        for (x) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;
        for (y) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

        sparse.spmv(mass, x, ax);
        sparse.spmv(mass, y, ay);

        var dot_axy: f64 = 0.0;
        var dot_xay: f64 = 0.0;
        var xtmx: f64 = 0.0;
        for (ax, y) |axi, yi| dot_axy += axi * yi;
        for (x, ay) |xi, ayi| dot_xay += xi * ayi;
        for (x, ax) |xi, axi| xtmx += xi * axi;

        try testing.expectApproxEqRel(dot_axy, dot_xay, 1e-12);
        try testing.expect(xtmx > 0.0);
    }
}
