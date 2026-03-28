//! Whitney 1-form mass matrix for the Galerkin Hodge star ★₁.
//!
//! The diagonal Hodge star ★₁ = diag(dual_length / length) is only a consistent
//! discretization when primal and dual edges are orthogonal. On meshes without
//! this orthogonality (e.g., barycentric dual), the diagonal ★₁ converges to
//! the wrong operator — producing ~8% systematic field error that doesn't
//! decrease with refinement.
//!
//! The Whitney/Galerkin mass matrix M₁ replaces the diagonal ★₁ with:
//!
//!   M₁(eᵢ, eⱼ) = ∫_Ω Wᵢ · Wⱼ dA
//!
//! where Wᵢ are Whitney 1-form basis functions. This is a sparse symmetric
//! positive-definite matrix that gives second-order field convergence on any
//! triangulation.
//!
//! Whitney 1-form for edge [a, b]:  W_{ab} = λₐ ∇λᵦ − λᵦ ∇λₐ
//!
//! The integral expands using the barycentric coordinate product rule:
//!   ∫_T λᵢ λⱼ dA = A(1 + δᵢⱼ) / 12
//!
//! and the fact that ∇λᵢ is constant on each triangle.
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

/// Assemble the Whitney 1-form mass matrix M₁ for a 2D simplicial mesh.
///
/// Returns a symmetric positive-definite CsrMatrix(f64) of size n_edges × n_edges.
/// Each triangle contributes a local 3×3 block to the three edges of its boundary.
///
/// The orientation signs from the boundary operator ∂₂ are accounted for:
/// the global Whitney form for edge eᵢ on triangle T is σᵢ · W_{local}, where
/// σᵢ = ±1 is the sign from ∂₂. The mass matrix entry gets multiplied by σᵢ · σⱼ.
pub fn assemble_whitney_mass_1(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const MeshType = @TypeOf(mesh.*);
    const n_edges = mesh.num_edges();
    const n_faces = mesh.num_faces();

    var assembler = sparse.TripletAssembler(f64).init(n_edges, n_edges);
    defer assembler.deinit(allocator);

    // Pre-fetch SoA slices for vertices, edges, and faces.
    const coords = mesh.vertices.slice().items(.coords);
    const mesh_edge_verts = mesh.edges.slice().items(.vertices);
    const face_verts = mesh.faces.slice().items(.vertices);

    // For each triangle, compute the local 3×3 Whitney mass matrix and
    // scatter into the global assembler.
    for (0..n_faces) |face_idx| {
        const fv = face_verts[face_idx];
        const v0 = coords[fv[0]];
        const v1 = coords[fv[1]];
        const v2 = coords[fv[2]];

        // Edge vectors opposite to each vertex.
        // e_opp[i] is the edge vector opposite to vertex i of the triangle.
        const e_opp = [3][MeshType.dimension]f64{
            vecSub(MeshType.dimension, v2, v1), // opposite v0: v1 → v2
            vecSub(MeshType.dimension, v0, v2), // opposite v1: v2 → v0
            vecSub(MeshType.dimension, v1, v0), // opposite v2: v0 → v1
        };

        // Triangle area (2D: half the cross product magnitude).
        const area = triangleArea(MeshType.dimension, v0, v1, v2);
        std.debug.assert(area > 0.0);

        // Gradient dot product matrix: G[i][j] = ∇λᵢ · ∇λⱼ = (e_opp[i] · e_opp[j]) / (4A²)
        const four_area_sq = 4.0 * area * area;
        var grad_dot: [3][3]f64 = undefined;
        for (0..3) |i| {
            for (0..3) |j| {
                grad_dot[i][j] = vecDot(MeshType.dimension, e_opp[i], e_opp[j]) / four_area_sq;
            }
        }

        // Local edges of the triangle in boundary orientation:
        // local edge 0: v0 → v1 (opposite v2)
        // local edge 1: v1 → v2 (opposite v0)
        // local edge 2: v2 → v0 (opposite v1)
        //
        // Each local edge [a, b] has Whitney form W_{ab} = λₐ∇λᵦ − λᵦ∇λₐ.
        // The vertex pairs (a, b) for each local edge:
        const local_edges = [3][2]u8{
            .{ 0, 1 }, // local edge 0: vertices 0 → 1
            .{ 1, 2 }, // local edge 1: vertices 1 → 2
            .{ 2, 0 }, // local edge 2: vertices 2 → 0
        };

        // Compute the local 3×3 mass matrix.
        // M_local(i, j) = ∫_T W_{aᵢbᵢ} · W_{aⱼbⱼ} dA
        //   = G[bᵢ][bⱼ]·S[aᵢ][aⱼ] − G[bᵢ][aⱼ]·S[aᵢ][bⱼ]
        //     − G[aᵢ][bⱼ]·S[bᵢ][aⱼ] + G[aᵢ][aⱼ]·S[bᵢ][bⱼ]
        // where S[i][j] = A·(1 + δᵢⱼ)/12
        var local_mass: [3][3]f64 = undefined;
        for (0..3) |i| {
            const ai = local_edges[i][0];
            const bi = local_edges[i][1];
            for (0..3) |j| {
                const aj = local_edges[j][0];
                const bj = local_edges[j][1];
                local_mass[i][j] = area * (grad_dot[bi][bj] * barycentricProduct(ai, aj) -
                    grad_dot[bi][aj] * barycentricProduct(ai, bj) -
                    grad_dot[ai][bj] * barycentricProduct(bi, aj) +
                    grad_dot[ai][aj] * barycentricProduct(bi, bj));
            }
        }

        // Map local edges to global edge indices and orientation signs
        // using the boundary operator ∂₂.
        const boundary_row = mesh.boundary_2.row(@intCast(face_idx));
        std.debug.assert(boundary_row.cols.len == 3);

        // boundary_2 row stores edges sorted by global index, which does NOT
        // in general match the local edge order (v0→v1, v1→v2, v2→v0).
        // Build an explicit mapping: local_edge_k → (global_edge, sign).
        // Match by looking up the global edge's vertex pair against the face's
        // vertex pair for each local edge.
        const boundary_global = boundary_row.cols;
        const boundary_signs = boundary_row.vals;

        var local_global_edge: [3]u32 = undefined;
        var local_global_sign: [3]i8 = undefined;
        for (0..3) |k| {
            const va = fv[local_edges[k][0]];
            const vb = fv[local_edges[k][1]];
            var found = false;
            for (0..3) |col| {
                const ev = mesh_edge_verts[boundary_global[col]];
                if ((ev[0] == va and ev[1] == vb) or (ev[0] == vb and ev[1] == va)) {
                    local_global_edge[k] = boundary_global[col];
                    local_global_sign[k] = boundary_signs[col];
                    found = true;
                    break;
                }
            }
            std.debug.assert(found); // every local edge must appear in ∂₂
        }

        // Scatter local mass into global assembler with orientation correction.
        for (0..3) |i| {
            for (0..3) |j| {
                const sign_i: f64 = @floatFromInt(local_global_sign[i]);
                const sign_j: f64 = @floatFromInt(local_global_sign[j]);
                const val = sign_i * sign_j * local_mass[i][j];
                try assembler.addEntry(allocator, local_global_edge[i], local_global_edge[j], val);
            }
        }
    }

    return assembler.build(allocator);
}

/// ∫_T λᵢ λⱼ dA / A = (1 + δᵢⱼ) / 12.
/// Returns the normalized barycentric product (divide by area separately).
fn barycentricProduct(i: u8, j: u8) f64 {
    return if (i == j) 1.0 / 6.0 else 1.0 / 12.0;
}

fn vecSub(comptime dim: usize, a: [dim]f64, b: [dim]f64) [dim]f64 {
    var result: [dim]f64 = undefined;
    for (0..dim) |i| result[i] = a[i] - b[i];
    return result;
}

fn vecDot(comptime dim: usize, a: [dim]f64, b: [dim]f64) f64 {
    var sum: f64 = 0;
    for (0..dim) |i| sum += a[i] * b[i];
    return sum;
}

fn triangleArea(comptime dim: usize, v0: [dim]f64, v1: [dim]f64, v2: [dim]f64) f64 {
    // 2D: A = ½|det([v1−v0, v2−v0])|
    if (dim == 2) {
        const dx1 = v1[0] - v0[0];
        const dy1 = v1[1] - v0[1];
        const dx2 = v2[0] - v0[0];
        const dy2 = v2[1] - v0[1];
        return @abs(dx1 * dy2 - dy1 * dx2) / 2.0;
    }
    // General dim: A = ½||(v1−v0) × (v2−v0)||
    @compileError("triangleArea not yet implemented for dim > 2");
}

const cg = @import("../math/cg.zig");

/// Precomputed Whitney Hodge star context for ★₁.
///
/// Holds the mass matrix, CG scratch space, and diagonal preconditioner.
/// Built once per mesh and reused across timesteps. The diagonal ★₁ values
/// (dual_length / length) serve as a spectrally equivalent preconditioner
/// for the CG solve, giving mesh-independent iteration counts.
pub fn WhitneyHodge1(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        mass_matrix: sparse.CsrMatrix(f64),
        scratch: cg.Scratch,
        precond_diagonal: []f64,

        /// Build the Whitney Hodge star context from a mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var mass = try assemble_whitney_mass_1(allocator, mesh);
            errdefer mass.deinit(allocator);

            var scratch = try cg.Scratch.init(allocator, mass.n_rows);
            errdefer scratch.deinit(allocator);

            // Build diagonal preconditioner from diagonal ★₁ = dual_length / length.
            const n_edges = mesh.num_edges();
            const precond = try allocator.alloc(f64, n_edges);
            errdefer allocator.free(precond);

            const edge_slice = mesh.edges.slice();
            const lengths = edge_slice.items(.length);
            const dual_lengths = edge_slice.items(.dual_length);
            for (precond, lengths, dual_lengths) |*p, len, dual_len| {
                p.* = dual_len / len;
            }

            return .{
                .mass_matrix = mass,
                .scratch = scratch,
                .precond_diagonal = precond,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.precond_diagonal);
            self.scratch.deinit(allocator);
            self.mass_matrix.deinit(allocator);
        }

        /// Apply ★₁: y = M₁ · x (forward SpMV).
        pub fn apply(self: Self, input: []const f64, output: []f64) void {
            sparse.spmv(self.mass_matrix, input, output);
        }

        /// Apply ★₁⁻¹: solve M₁ · x = b via preconditioned CG.
        /// Returns the solution in `x`. Uses `b` as the RHS.
        pub fn solve_inverse(self: *Self, b: []const f64, x: []f64) void {
            var precond = cg.DiagonalPreconditioner{ .diagonal = self.precond_diagonal };
            const result = cg.solve(
                self.mass_matrix,
                b,
                x,
                1e-10,
                1000,
                cg.DiagonalPreconditioner.apply,
                @ptrCast(&precond),
                self.scratch,
            );
            std.debug.assert(result.converged);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);

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
