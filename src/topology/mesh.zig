//! Simplicial mesh topology and geometry.
//!
//! Provides `Mesh(embedding_dimension, topological_dimension)`, a simplicial
//! complex with SoA entity storage, oriented boundary operators ∂ₖ for all
//! 1 ≤ k ≤ topological_dimension in CSR format, and geometric data (lengths,
//! areas, dual volumes).

const std = @import("std");
const testing = std.testing;
const flux = @import("../root.zig");
const sparse = @import("../math/sparse.zig");
const whitney = @import("../operators/whitney_mass.zig");

/// Boundary operators use packed sign storage: one bit per stored incidence.
pub const BoundaryMatrix = sparse.PackedIncidenceMatrix;
const DenseBoundaryMatrix = sparse.CsrMatrix(i8);

/// Standalone 0-simplex (vertex) type constructor.
///
/// Carries embedding coordinates and the volume of the dual cell
/// (Voronoi/circumcentric region).
pub fn Vertex(comptime mesh_embedding_dimension: usize) type {
    return struct {
        coords: [mesh_embedding_dimension]f64,
        /// Volume of the dual cell (cotangent-weighted).
        /// For topological_dimension=2 this is an area; for
        /// topological_dimension=3 a volume. The name is
        /// dimension-agnostic.
        dual_volume: f64,
    };
}

/// Standalone k-simplex type constructor.
///
/// Used both by mesh storage internals and by external code that needs the
/// concrete simplex record type for a given embedding/topological dimension
/// and simplex degree `k`.
///
/// A k-simplex (k ≥ 1) is defined by `k + 1` oriented vertex indices and
/// carries its k-dimensional measure and barycenter.
pub fn Simplex(comptime embedding_dimension: usize, comptime topological_dimension: usize, comptime k: usize) type {
    if (k < 1) @compileError("Simplex(k) requires k ≥ 1; vertices are stored separately");
    if (k > topological_dimension) @compileError(std.fmt.comptimePrint(
        "Simplex({d}) exceeds topological dimension {d}",
        .{ k, topological_dimension },
    ));
    return struct {
        /// Vertex indices defining this oriented k-simplex.
        vertices: [k + 1]u32,
        /// k-dimensional measure (length for k=1, area for k=2, volume for k=3).
        volume: f64,
        /// Barycenter (centroid) coordinates.
        barycenter: [embedding_dimension]f64,
    };
}

// ───────────────────────────────────────────────────────────────────────────
// Comptime simplex storage type
// ───────────────────────────────────────────────────────────────────────────

/// Build a tuple type whose field `i` is
/// `MultiArrayList(Simplex(embedding_dimension, topological_dimension, i + 1))`
/// for `i = 0..topological_dimension-1`.
///
/// This lets us store all k-simplex lists (`k = 1..topological_dimension`) in
/// a single comptime-indexed struct without hardcoded field names.
fn SimplexListsType(comptime embedding_dimension: usize, comptime topological_dimension: usize) type {
    var fields: [topological_dimension]std.builtin.Type.StructField = undefined;
    for (0..topological_dimension) |i| {
        const k = i + 1;
        const MAL = std.MultiArrayList(Simplex(embedding_dimension, topological_dimension, k));
        fields[i] = .{
            .name = std.fmt.comptimePrint("{d}", .{i}),
            .type = MAL,
            .default_value_ptr = null,
            .is_comptime = false,
            .alignment = @alignOf(MAL),
        };
    }
    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &fields,
        .decls = &.{},
        .is_tuple = true,
    } });
}

// ───────────────────────────────────────────────────────────────────────────
// Mesh
// ───────────────────────────────────────────────────────────────────────────

/// Simplicial mesh parameterized on embedding and topological dimension.
///
/// `embedding_dimension`: the ambient space ℝⁿ that vertices live in.
/// `topological_dimension`: the intrinsic dimension of the simplicial complex
/// (2 for triangulated surfaces, 3 for tetrahedral volumes).
///
/// Topological entities use standalone `Vertex(embedding_dimension)` and
/// `Simplex(embedding_dimension, topological_dimension, k)` record types,
/// stored in `std.MultiArrayList` for SoA cache-friendly layout.
///
/// Boundary operators ∂ₖ for 1 ≤ k ≤ topological_dimension are stored in CSR
/// format.
pub fn Mesh(comptime mesh_embedding_dimension: usize, comptime mesh_topological_dimension: usize) type {
    comptime {
        if (mesh_embedding_dimension < 1) @compileError("embedding dimension must be at least 1");
        if (mesh_topological_dimension < 2) @compileError("topological dimension must be at least 2");
        if (mesh_topological_dimension > 3) @compileError("topological dimension > 3 is not yet supported");
        if (mesh_topological_dimension > mesh_embedding_dimension) @compileError("topological dimension cannot exceed embedding dimension");
    }

    return struct {
        const Self = @This();

        /// Embedding dimension — the ambient space ℝⁿ.
        pub const embedding_dimension = mesh_embedding_dimension;

        /// Topological dimension of the simplicial complex.
        pub const topological_dimension = mesh_topological_dimension;

        /// Comptime-generated tuple of
        /// `MultiArrayList(Simplex(mesh_embedding_dimension, mesh_topological_dimension, k))`
        /// for `k = 1..topological_dimension`.
        /// Access via `simplices(k)` which returns a SoA slice.
        const SimplexLists = SimplexListsType(mesh_embedding_dimension, mesh_topological_dimension);

        // -- Storage --

        vertices: std.MultiArrayList(Vertex(mesh_embedding_dimension)),

        /// SoA storage for k-simplices (k = 1..topological_dimension).
        /// Indexed as `simplex_lists[k-1]`.
        /// Use `simplices(k)` for typed access.
        simplex_lists: SimplexLists,

        /// Boundary operators ∂ₖ for k = 1..topological_dimension, stored as
        /// `boundaries[k-1]`.
        /// ∂₁: edges → vertices, ∂₂: faces → edges, ∂₃: tets → faces.
        boundaries: [topological_dimension]BoundaryMatrix,

        /// Dual cell volumes for each 1-simplex (edge). For
        /// topological_dimension=2 this is the dual edge length; for
        /// topological_dimension=3 the dual face area. Stored separately
        /// from Simplex(1) to avoid every simplex carrying dimension-dependent fields.
        dual_edge_volumes: []f64,

        /// Indices of edges on the mesh boundary (adjacent to exactly one face).
        /// Precomputed during construction for use by boundary condition routines.
        boundary_edges: []u32,

        /// Whitney 1-form mass matrix M₁ for the Galerkin Hodge star ★₁.
        /// M₁(eᵢ, eⱼ) = ∫_Ω Wᵢ · Wⱼ dA — sparse SPD, exact on any triangulation.
        /// Used by hodge_star for k=1 (SpMV) and hodge_star_inverse for k=1 (CG solve).
        whitney_mass_1: sparse.CsrMatrix(f64),

        /// Diagonal preconditioner for CG solves with whitney_mass_1.
        /// Values are the diagonal ★₁ = dual_volume / volume per edge —
        /// spectrally equivalent to M₁, giving mesh-independent iteration counts.
        preconditioner_1: []f64,

        // -- Lifetime --

        /// Free all entity storage and boundary matrices.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.preconditioner_1);
            self.whitney_mass_1.deinit(allocator);
            allocator.free(self.boundary_edges);
            allocator.free(self.dual_edge_volumes);
            self.vertices.deinit(allocator);
            inline for (0..topological_dimension) |i| {
                self.simplex_lists[i].deinit(allocator);
            }
            for (&self.boundaries) |*b| b.deinit(allocator);
        }

        // -- Accessors --

        /// Return the boundary operator ∂ₖ for the given degree
        /// (1 ≤ k ≤ topological_dimension).
        ///
        /// ∂₁ maps edges → vertices, ∂₂ maps faces → edges, ∂₃ maps tets → faces.
        /// Stored in coboundary orientation (rows indexed by higher-dimensional
        /// cells), so ∂ₖ also serves as the exterior derivative dₖ₋₁.
        pub fn boundary(self: Self, comptime k: comptime_int) BoundaryMatrix {
            if (k < 1 or k > topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no boundary operator ∂_{d} on a {d}-dimensional mesh (need 1 ≤ k ≤ {d})",
                    .{ k, topological_dimension, topological_dimension },
                ));
            }
            return self.boundaries[k - 1];
        }

        /// SoA slice of k-simplices (`k ≥ 1`).
        ///
        /// Returns a
        /// `std.MultiArrayList(Simplex(embedding_dimension, topological_dimension, k)).Slice`
        /// with fields `.vertices`, `.volume`, and `.barycenter`.
        pub fn simplices(self: *const Self, comptime k: comptime_int) std.MultiArrayList(Simplex(embedding_dimension, topological_dimension, k)).Slice {
            if (k < 1 or k > topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no {d}-simplices on a {d}-dimensional mesh (need 1 ≤ k ≤ {d})",
                    .{ k, topological_dimension, topological_dimension },
                ));
            }
            return self.simplex_lists[k - 1].slice();
        }

        /// Number of k-simplices in the mesh (generic entity count accessor).
        pub fn num_cells(self: Self, comptime k: comptime_int) u32 {
            if (k == 0) return @intCast(self.vertices.len);
            if (k < 0 or k > topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no {d}-cells on a {d}-dimensional mesh",
                    .{ k, topological_dimension },
                ));
            }
            return @intCast(self.simplex_lists[k - 1].len);
        }

        /// Number of vertices in the mesh.
        pub fn num_vertices(self: Self) u32 {
            return self.num_cells(0);
        }

        /// Number of edges in the mesh.
        pub fn num_edges(self: Self) u32 {
            return self.num_cells(1);
        }

        /// Number of faces in the mesh.
        pub fn num_faces(self: Self) u32 {
            return self.num_cells(2);
        }

        /// Number of tetrahedra in the mesh (only for topological_dimension ≥ 3).
        pub fn num_tets(self: Self) u32 {
            return self.num_cells(3);
        }

        // ───────────────────────────────────────────────────────────────────
        // Uniform grid constructor
        // ───────────────────────────────────────────────────────────────────

        /// Construct a uniform triangulated rectangular grid on `[0, width] × [0, height]`.
        ///
        /// Only available for 2D meshes (`topological_dimension = 2`). For
        /// 3D meshes, use the
        /// tetrahedral grid constructor (see #81).
        ///
        /// The rectangle is divided into `nx × ny` rectangular cells, each split
        /// by a SW-to-NE diagonal into two CCW-oriented triangles (`2 * nx * ny` total).
        /// All geometric data — lengths, areas, barycenters, dual volumes — is
        /// computed during construction.
        ///
        /// Uses a **barycentric dual** mesh: dual vertices are placed at triangle
        /// centroids rather than circumcenters. This avoids the degeneracy where
        /// right triangles sharing a hypotenuse have coincident circumcenters,
        /// ensuring every edge has nonzero dual length and ★₁ is non-singular.
        pub fn uniform_grid(
            allocator: std.mem.Allocator,
            nx: u32,
            ny: u32,
            width: f64,
            height: f64,
        ) !Self {
            comptime {
                if (topological_dimension != 2) @compileError("uniform_grid is only available for 2D meshes (topological_dimension = 2)");
            }
            if (nx == 0 or ny == 0) @panic("grid dimensions must be positive");
            if (!(width > 0) or !(height > 0)) @panic("domain size must be positive");

            const dx = width / @as(f64, @floatFromInt(nx));
            const dy = height / @as(f64, @floatFromInt(ny));

            const vertex_count: u32 = (nx + 1) * (ny + 1);
            const horizontal_edge_count: u32 = nx * (ny + 1);
            const vertical_edge_count: u32 = (nx + 1) * ny;
            const diagonal_edge_count: u32 = nx * ny;
            const edge_count: u32 = horizontal_edge_count + vertical_edge_count + diagonal_edge_count;
            const face_count: u32 = 2 * nx * ny;

            // -- Allocate entity storage --

            var vertices = std.MultiArrayList(Vertex(embedding_dimension)){};
            try vertices.ensureTotalCapacity(allocator, vertex_count);
            errdefer vertices.deinit(allocator);

            var edges_list = std.MultiArrayList(Simplex(embedding_dimension, topological_dimension, 1)){};
            try edges_list.ensureTotalCapacity(allocator, edge_count);
            errdefer edges_list.deinit(allocator);

            var faces_list = std.MultiArrayList(Simplex(embedding_dimension, topological_dimension, 2)){};
            try faces_list.ensureTotalCapacity(allocator, face_count);
            errdefer faces_list.deinit(allocator);

            // -- Populate vertices --
            // Indexed as vertex_index(i, j) = i * (ny + 1) + j.
            for (0..nx + 1) |i_u| {
                const fi: f64 = @floatFromInt(i_u);
                for (0..ny + 1) |j_u| {
                    const fj: f64 = @floatFromInt(j_u);
                    var coords: [embedding_dimension]f64 = @splat(0);
                    coords[0] = fi * dx;
                    coords[1] = fj * dy;
                    vertices.appendAssumeCapacity(.{ .coords = coords, .dual_volume = 0 });
                }
            }

            // -- Populate edges --
            // Horizontal: horizontal_edge(i,j) = j*nx + i, tail = vertex(i,j), head = vertex(i+1,j)
            for (0..ny + 1) |j_u| {
                for (0..nx) |i_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i + 1, j, ny) },
                        .volume = 0,
                        .barycenter = @splat(0),
                    });
                }
            }
            // Vertical: vertical_edge(i,j) = horizontal_edge_count + i*ny + j
            for (0..nx + 1) |i_u| {
                for (0..ny) |j_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i, j + 1, ny) },
                        .volume = 0,
                        .barycenter = @splat(0),
                    });
                }
            }
            // Diagonal (SW→NE): diagonal_edge(i,j) = horizontal_edge_count + vertical_edge_count + i*ny + j
            for (0..nx) |i_u| {
                for (0..ny) |j_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i + 1, j + 1, ny) },
                        .volume = 0,
                        .barycenter = @splat(0),
                    });
                }
            }

            // -- Populate faces --
            // Two triangles per cell, both CCW-oriented.
            for (0..nx) |i_u| {
                const i: u32 = @intCast(i_u);
                for (0..ny) |j_u| {
                    const j: u32 = @intCast(j_u);
                    const sw = vertex_index(i, j, ny);
                    const se = vertex_index(i + 1, j, ny);
                    const nw = vertex_index(i, j + 1, ny);
                    const ne = vertex_index(i + 1, j + 1, ny);

                    const zero_cc: [embedding_dimension]f64 = @splat(0);

                    // Lower-right triangle: SW → SE → NE
                    faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, se, ne }, .volume = 0, .barycenter = zero_cc });
                    // Upper-left triangle:  SW → NE → NW
                    faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, ne, nw }, .volume = 0, .barycenter = zero_cc });
                }
            }

            // -- Build ∂₁ (n_edges × n_vertices) --
            // Each edge row has exactly 2 nonzeros: tail = −1, head = +1.
            var boundary_1: BoundaryMatrix = undefined;
            {
                var boundary_1_dense = try DenseBoundaryMatrix.init(allocator, edge_count, vertex_count, 2 * edge_count);
                errdefer boundary_1_dense.deinit(allocator);
                const edge_verts = edges_list.slice().items(.vertices);
                for (0..edge_count) |e| {
                    boundary_1_dense.row_ptr[e] = @intCast(2 * e);
                    // tail < head for the uniform grid, so columns are sorted
                    boundary_1_dense.col_idx[2 * e] = edge_verts[e][0];
                    boundary_1_dense.values[2 * e] = -1;
                    boundary_1_dense.col_idx[2 * e + 1] = edge_verts[e][1];
                    boundary_1_dense.values[2 * e + 1] = 1;
                }
                boundary_1_dense.row_ptr[edge_count] = 2 * edge_count;

                boundary_1 = try BoundaryMatrix.fromBoundaryCsr(allocator, 1, boundary_1_dense);
                boundary_1_dense.deinit(allocator);
            }
            errdefer boundary_1.deinit(allocator);

            // -- Build ∂₂ (n_faces × n_edges) --
            // Each face row has exactly 3 nonzeros for the oriented boundary edges.
            //
            // Lower-right (SW,SE,NE): +h(i,j), +vert(i+1,j), −diag(i,j)
            // Upper-left  (SW,NE,NW): −h(i,j+1), −vert(i,j), +diag(i,j)
            //
            // Column ordering: horizontal < vertical < diagonal (always sorted).
            var boundary_2: BoundaryMatrix = undefined;
            {
                var boundary_2_dense = try DenseBoundaryMatrix.init(allocator, face_count, edge_count, 3 * face_count);
                errdefer boundary_2_dense.deinit(allocator);
                var f_idx: u32 = 0;
                for (0..nx) |i_u| {
                    const i: u32 = @intCast(i_u);
                    for (0..ny) |j_u| {
                        const j: u32 = @intCast(j_u);
                        const h_ij = horizontal_edge_index(i, j, nx);
                        const h_i_jp1 = horizontal_edge_index(i, j + 1, nx);
                        const v_ip1_j = vertical_edge_index(i + 1, j, ny, horizontal_edge_count);
                        const v_i_j = vertical_edge_index(i, j, ny, horizontal_edge_count);
                        const d_ij = diagonal_edge_index(i, j, ny, horizontal_edge_count, vertical_edge_count);

                        // Lower-right triangle
                        boundary_2_dense.row_ptr[f_idx] = 3 * f_idx;
                        boundary_2_dense.col_idx[3 * f_idx + 0] = h_ij;
                        boundary_2_dense.values[3 * f_idx + 0] = 1;
                        boundary_2_dense.col_idx[3 * f_idx + 1] = v_ip1_j;
                        boundary_2_dense.values[3 * f_idx + 1] = 1;
                        boundary_2_dense.col_idx[3 * f_idx + 2] = d_ij;
                        boundary_2_dense.values[3 * f_idx + 2] = -1;
                        f_idx += 1;

                        // Upper-left triangle
                        boundary_2_dense.row_ptr[f_idx] = 3 * f_idx;
                        boundary_2_dense.col_idx[3 * f_idx + 0] = h_i_jp1;
                        boundary_2_dense.values[3 * f_idx + 0] = -1;
                        boundary_2_dense.col_idx[3 * f_idx + 1] = v_i_j;
                        boundary_2_dense.values[3 * f_idx + 1] = -1;
                        boundary_2_dense.col_idx[3 * f_idx + 2] = d_ij;
                        boundary_2_dense.values[3 * f_idx + 2] = 1;
                        f_idx += 1;
                    }
                }
                boundary_2_dense.row_ptr[face_count] = 3 * face_count;

                boundary_2 = try BoundaryMatrix.fromBoundaryCsr(allocator, 2, boundary_2_dense);
                boundary_2_dense.deinit(allocator);
            }
            errdefer boundary_2.deinit(allocator);

            // -- Compute primal geometry --

            const coords = vertices.slice().items(.coords);

            // Edge lengths and barycenters
            {
                const edge_slice = edges_list.slice();
                const edge_verts = edge_slice.items(.vertices);
                const edge_volumes = edge_slice.items(.volume);
                const edge_barycenters = edge_slice.items(.barycenter);
                for (0..edge_count) |e| {
                    edge_volumes[e] = euclidean_distance(coords[edge_verts[e][0]], coords[edge_verts[e][1]]);
                    edge_barycenters[e] = point_midpoint(coords[edge_verts[e][0]], coords[edge_verts[e][1]]);
                }
            }

            // Face areas and barycenters
            {
                const face_slice = faces_list.slice();
                const face_verts = face_slice.items(.vertices);
                const face_volumes = face_slice.items(.volume);
                const face_barycenters = face_slice.items(.barycenter);
                for (0..face_count) |f| {
                    const p0 = coords[face_verts[f][0]];
                    const p1 = coords[face_verts[f][1]];
                    const p2 = coords[face_verts[f][2]];
                    face_volumes[f] = triangle_area(p0, p1, p2);
                    face_barycenters[f] = triangle_barycenter(p0, p1, p2);
                }
            }

            // -- Compute circumcentric dual geometry --

            // Dual edge volumes require edge→face adjacency.
            // Each edge borders at most 2 faces; boundary edges border exactly 1.
            // Single scratch allocation: [count | face_0 | face_1], each edge_count u32s.
            const dual_edge_volumes = try allocator.alloc(f64, edge_count);
            errdefer allocator.free(dual_edge_volumes);
            var boundary_edge_buf: []u32 = &.{};
            {
                const scratch = try allocator.alloc(u32, 3 * edge_count);
                defer allocator.free(scratch);
                const edge_face_count = scratch[0..edge_count];
                const edge_face_0 = scratch[edge_count .. 2 * edge_count];
                const edge_face_1 = scratch[2 * edge_count .. 3 * edge_count];
                @memset(edge_face_count, 0);

                for (0..face_count) |f| {
                    const face_edges = boundary_2.row(@intCast(f));
                    for (face_edges.cols) |e| {
                        const count = edge_face_count[e];
                        if (count == 0) {
                            edge_face_0[e] = @intCast(f);
                        } else {
                            edge_face_1[e] = @intCast(f);
                        }
                        edge_face_count[e] = count + 1;
                    }
                }

                const barycenters = faces_list.slice().items(.barycenter);
                const edge_verts = edges_list.slice().items(.vertices);

                // Count boundary edges (adjacent to exactly one face) and
                // compute dual edge volumes in one pass.
                var boundary_count: u32 = 0;
                for (0..edge_count) |e| {
                    if (edge_face_count[e] == 2) {
                        dual_edge_volumes[e] = euclidean_distance(
                            barycenters[edge_face_0[e]],
                            barycenters[edge_face_1[e]],
                        );
                    } else if (edge_face_count[e] == 1) {
                        const mid = point_midpoint(coords[edge_verts[e][0]], coords[edge_verts[e][1]]);
                        dual_edge_volumes[e] = euclidean_distance(barycenters[edge_face_0[e]], mid);
                        boundary_count += 1;
                    } else {
                        return flux.Error.NonManifoldEdge;
                    }
                }

                // Collect boundary edge indices.
                boundary_edge_buf = try allocator.alloc(u32, boundary_count);
                var bi: u32 = 0;
                for (0..edge_count) |e| {
                    if (edge_face_count[e] == 1) {
                        boundary_edge_buf[bi] = @intCast(e);
                        bi += 1;
                    }
                }
                std.debug.assert(bi == boundary_count);
            }

            // Dual vertex volumes via the cotangent formula.
            //
            // For triangle (v₀, v₁, v₂) with edge lengths l₀₁, l₁₂, l₂₀:
            //   cot(αᵢ) = (lⱼₖ² + lₖᵢ² − lᵢⱼ²) / (4 · area)
            //
            //   dual_volume[vᵢ] += (cot(αⱼ) · lᵢₖ² + cot(αₖ) · lᵢⱼ²) / 8
            //
            // This gives the circumcentric (Voronoi) dual area contribution.
            // Reference: Meyer et al., "Discrete Differential-Geometry Operators
            // for Triangulated 2-Manifolds" (2002), §4.
            {
                const dual_volumes = vertices.slice().items(.dual_volume);
                @memset(dual_volumes, 0);

                const face_verts = faces_list.slice().items(.vertices);
                for (0..face_count) |f| {
                    const vs = face_verts[f];
                    const p0 = coords[vs[0]];
                    const p1 = coords[vs[1]];
                    const p2 = coords[vs[2]];

                    const l01_sq = distance_squared(p0, p1);
                    const l12_sq = distance_squared(p1, p2);
                    const l20_sq = distance_squared(p2, p0);

                    const area_4 = 4.0 * triangle_area(p0, p1, p2);
                    // Degenerate triangles cause division by zero in the
                    // cotangent formula below.
                    if (!(area_4 > 0)) return flux.Error.DegenerateTriangle;

                    const cot0 = (l01_sq + l20_sq - l12_sq) / area_4;
                    const cot1 = (l01_sq + l12_sq - l20_sq) / area_4;
                    const cot2 = (l12_sq + l20_sq - l01_sq) / area_4;

                    dual_volumes[vs[0]] += (cot1 * l20_sq + cot2 * l01_sq) / 8.0;
                    dual_volumes[vs[1]] += (cot2 * l01_sq + cot0 * l12_sq) / 8.0;
                    dual_volumes[vs[2]] += (cot0 * l12_sq + cot1 * l20_sq) / 8.0;
                }
            }

            // -- Whitney mass matrix M₁ and diagonal preconditioner --
            //
            // M₁ is the Galerkin inner product of Whitney 1-form basis functions.
            // It replaces the diagonal ★₁ (which only converges on orthogonal
            // dual meshes) with an SPD matrix that is exact on any triangulation.
            // The diagonal ★₁ values (dual_length / length) serve as a spectrally
            // equivalent preconditioner for CG solves of M₁.

            var partial = Self{
                .vertices = vertices,
                .simplex_lists = .{ edges_list, faces_list },
                .boundaries = .{ boundary_1, boundary_2 },
                .dual_edge_volumes = dual_edge_volumes,
                .boundary_edges = boundary_edge_buf,
                .whitney_mass_1 = undefined,
                .preconditioner_1 = undefined,
            };

            var mass = try whitney.assemble_whitney_mass_1(allocator, &partial);
            errdefer mass.deinit(allocator);

            const precond = try allocator.alloc(f64, edge_count);
            errdefer allocator.free(precond);
            {
                const edge_volumes = edges_list.slice().items(.volume);
                for (precond, edge_volumes, dual_edge_volumes) |*p, volume, dual_volume| {
                    p.* = dual_volume / volume;
                }
            }

            partial.whitney_mass_1 = mass;
            partial.preconditioner_1 = precond;
            return partial;
        }

        // ───────────────────────────────────────────────────────────────────
        // Index helpers (grid-specific)
        // ───────────────────────────────────────────────────────────────────

        fn vertex_index(i: u32, j: u32, ny_val: u32) u32 {
            return i * (ny_val + 1) + j;
        }

        fn horizontal_edge_index(i: u32, j: u32, nx_val: u32) u32 {
            return j * nx_val + i;
        }

        fn vertical_edge_index(i: u32, j: u32, ny_val: u32, n_horizontal: u32) u32 {
            return n_horizontal + i * ny_val + j;
        }

        fn diagonal_edge_index(i: u32, j: u32, ny_val: u32, n_horizontal: u32, n_vertical: u32) u32 {
            return n_horizontal + n_vertical + i * ny_val + j;
        }

        // ───────────────────────────────────────────────────────────────────
        // Geometry primitives
        // ───────────────────────────────────────────────────────────────────

        fn distance_squared(a: [embedding_dimension]f64, b: [embedding_dimension]f64) f64 {
            var sum: f64 = 0;
            inline for (0..embedding_dimension) |d| {
                const diff = b[d] - a[d];
                sum += diff * diff;
            }
            return sum;
        }

        fn euclidean_distance(a: [embedding_dimension]f64, b: [embedding_dimension]f64) f64 {
            return @sqrt(distance_squared(a, b));
        }

        fn point_midpoint(a: [embedding_dimension]f64, b: [embedding_dimension]f64) [embedding_dimension]f64 {
            var result: [embedding_dimension]f64 = undefined;
            inline for (0..embedding_dimension) |d| {
                result[d] = 0.5 * (a[d] + b[d]);
            }
            return result;
        }

        /// Signed area of a 2D triangle (positive for CCW orientation).
        /// Uses only the first two coordinate components.
        fn signed_triangle_area(p0: [embedding_dimension]f64, p1: [embedding_dimension]f64, p2: [embedding_dimension]f64) f64 {
            const dx1 = p1[0] - p0[0];
            const dy1 = p1[1] - p0[1];
            const dx2 = p2[0] - p0[0];
            const dy2 = p2[1] - p0[1];
            return (dx1 * dy2 - dx2 * dy1) / 2.0;
        }

        /// Unsigned area of a 2D triangle. Uses only the first two coordinate components.
        fn triangle_area(p0: [embedding_dimension]f64, p1: [embedding_dimension]f64, p2: [embedding_dimension]f64) f64 {
            return @abs(signed_triangle_area(p0, p1, p2));
        }

        /// Barycenter (centroid) of a triangle: average of its three vertices.
        fn triangle_barycenter(a: [embedding_dimension]f64, b: [embedding_dimension]f64, c: [embedding_dimension]f64) [embedding_dimension]f64 {
            var result: [embedding_dimension]f64 = undefined;
            inline for (0..embedding_dimension) |d| {
                result[d] = (a[d] + b[d] + c[d]) / 3.0;
            }
            return result;
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

test "uniform grid entity counts" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(u32, 20), mesh.num_vertices()); // (3+1)*(4+1)
    try testing.expectEqual(@as(u32, 43), mesh.num_edges()); // 3*5 + 4*4 + 3*4
    try testing.expectEqual(@as(u32, 24), mesh.num_faces()); // 2*3*4
}

test "boundary operator ∂₁ has exactly 2 nonzeros per row" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    for (0..mesh.num_edges()) |e| {
        const r = mesh.boundary(1).row(@intCast(e));
        try testing.expectEqual(@as(usize, 2), r.cols.len);
        // tail = −1, head = +1
        try testing.expectEqual(@as(i8, -1), r.sign(0));
        try testing.expectEqual(@as(i8, 1), r.sign(1));
        // columns sorted
        try testing.expect(r.cols[0] < r.cols[1]);
    }
}

test "boundary operator ∂₂ has exactly 3 nonzeros per row" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    for (0..mesh.num_faces()) |f| {
        const r = mesh.boundary(2).row(@intCast(f));
        try testing.expectEqual(@as(usize, 3), r.cols.len);
        // columns sorted
        try testing.expect(r.cols[0] < r.cols[1]);
        try testing.expect(r.cols[1] < r.cols[2]);
    }
}

test "boundary of boundary is zero for 2D triangulations" {
    // ∂₁ ∘ ∂₂ = 0: for each face, computing the boundary of its boundary
    // at every vertex must yield zero. This is tested with integer arithmetic
    // (exact, no floating-point tolerance).
    const allocator = testing.allocator;

    const sizes = [_][2]u32{
        .{ 1, 1 },   .{ 1, 2 },  .{ 2, 1 },  .{ 2, 2 },
        .{ 3, 3 },   .{ 4, 5 },  .{ 5, 4 },  .{ 7, 3 },
        .{ 10, 10 }, .{ 1, 10 }, .{ 10, 1 },
    };

    for (sizes) |size| {
        var mesh = try Mesh(2, 2).uniform_grid(allocator, size[0], size[1], 1.0, 1.0);
        defer mesh.deinit(allocator);

        var vertex_sum = try allocator.alloc(i32, mesh.num_vertices());
        defer allocator.free(vertex_sum);

        for (0..mesh.num_faces()) |f| {
            @memset(vertex_sum, 0);

            const face_row = mesh.boundary(2).row(@intCast(f));
            for (face_row.cols, 0..) |edge_idx, face_entry_idx| {
                const face_sign = face_row.sign(face_entry_idx);
                const edge_row = mesh.boundary(1).row(edge_idx);
                for (edge_row.cols, 0..) |vert_idx, edge_entry_idx| {
                    const edge_sign = edge_row.sign(edge_entry_idx);
                    vertex_sum[vert_idx] += @as(i32, face_sign) * @as(i32, edge_sign);
                }
            }

            for (vertex_sum) |s| {
                try testing.expectEqual(@as(i32, 0), s);
            }
        }
    }
}

test "edge lengths for unit square 1×1 grid" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const lengths = mesh.simplices(1).items(.volume);

    // 1×1 grid: 2 horizontal (len 1), 2 vertical (len 1), 1 diagonal (len √2)
    try testing.expectEqual(@as(u32, 5), mesh.num_edges());
    // Horizontal edges
    try testing.expectApproxEqAbs(@as(f64, 1.0), lengths[0], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 1.0), lengths[1], 1e-15);
    // Vertical edges
    try testing.expectApproxEqAbs(@as(f64, 1.0), lengths[2], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 1.0), lengths[3], 1e-15);
    // Diagonal edge
    try testing.expectApproxEqAbs(@sqrt(2.0), lengths[4], 1e-15);
}

test "face areas for uniform grid" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const areas = mesh.simplices(2).items(.volume);
    const expected_area: f64 = 0.5 * 0.5 / 2.0; // dx * dy / 2

    for (0..mesh.num_faces()) |f| {
        try testing.expectApproxEqAbs(expected_area, areas[f], 1e-15);
    }
}

test "barycenters of triangles are at centroids" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 1, 1, 2.0, 2.0);
    defer mesh.deinit(allocator);

    const barycenters = mesh.simplices(2).items(.barycenter);

    // Single 2×2 cell split by SW→NE diagonal:
    //   lower triangle: (0,0), (2,0), (2,2) → barycenter (4/3, 2/3)
    //   upper triangle: (0,0), (2,2), (0,2) → barycenter (2/3, 4/3)
    try testing.expectApproxEqAbs(4.0 / 3.0, barycenters[0][0], 1e-15);
    try testing.expectApproxEqAbs(2.0 / 3.0, barycenters[0][1], 1e-15);
    try testing.expectApproxEqAbs(2.0 / 3.0, barycenters[1][0], 1e-15);
    try testing.expectApproxEqAbs(4.0 / 3.0, barycenters[1][1], 1e-15);
}

test "dual vertex areas sum to total mesh area" {
    const allocator = testing.allocator;
    const width: f64 = 3.0;
    const height: f64 = 2.0;

    var mesh = try Mesh(2, 2).uniform_grid(allocator, 5, 4, width, height);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_volume);
    var total: f64 = 0;
    for (dual_areas) |a| {
        total += a;
    }

    try testing.expectApproxEqAbs(width * height, total, 1e-12);
}

test "interior vertex dual area equals cell area" {
    const allocator = testing.allocator;
    const nx: u32 = 5;
    const ny: u32 = 4;
    const width: f64 = 3.0;
    const height: f64 = 2.0;
    const dx = width / @as(f64, @floatFromInt(nx));
    const dy = height / @as(f64, @floatFromInt(ny));

    var mesh = try Mesh(2, 2).uniform_grid(allocator, nx, ny, width, height);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_volume);

    // Interior vertices should each have dual area = dx * dy
    for (1..nx) |i| {
        for (1..ny) |j| {
            const idx = Mesh(2, 2).vertex_index(@intCast(i), @intCast(j), ny);
            try testing.expectApproxEqAbs(dx * dy, dual_areas[idx], 1e-13);
        }
    }
}

test "all edges have nonzero dual length" {
    // With barycentric dual, every edge (including diagonals) has nonzero
    // dual length because adjacent triangles have distinct barycenters.
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const dual_lengths = mesh.dual_edge_volumes;

    for (0..mesh.num_edges()) |e| {
        try testing.expect(dual_lengths[e] > 0.0);
    }
}

test "Mesh(2, 2) compiles at dimension 2" {
    // Compile-time check: Mesh(2, 2) is a valid type.
    const M = Mesh(2, 2);
    try testing.expect(M.embedding_dimension == 2);
}

test "Mesh(3, 2) compiles — surface in ℝ³" {
    const M = Mesh(3, 2);
    try testing.expect(M.embedding_dimension == 3);
    try testing.expect(M.topological_dimension == 2);
}

test "Mesh(3, 3) compiles — volume in ℝ³" {
    const M = Mesh(3, 3);
    try testing.expect(M.embedding_dimension == 3);
    try testing.expect(M.topological_dimension == 3);
}

test "random grid dimensions produce valid meshes (100 trials)" {
    // Stress test: random (nx, ny, width, height) tuples must all produce
    // meshes with valid geometry — positive edge lengths, positive face areas,
    // positive dual areas that sum to the total domain area, and nonzero
    // dual lengths on all edges.
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0xE5B_57E55);

    for (0..100) |_| {
        const nx: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 20));
        const ny: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 20));
        const width = rng.random().float(f64) * 99.9 + 0.1; // (0.1, 100.0)
        const height = rng.random().float(f64) * 99.9 + 0.1;

        var mesh = try Mesh(2, 2).uniform_grid(allocator, nx, ny, width, height);
        defer mesh.deinit(allocator);

        // Entity counts match the grid formulas.
        try testing.expectEqual((nx + 1) * (ny + 1), mesh.num_vertices());
        try testing.expectEqual(2 * nx * ny, mesh.num_faces());

        // All edge lengths are positive.
        const lengths = mesh.simplices(1).items(.volume);
        for (lengths) |l| {
            try testing.expect(l > 0.0);
        }

        // All face areas are positive.
        const areas = mesh.simplices(2).items(.volume);
        for (areas) |a| {
            try testing.expect(a > 0.0);
        }

        // All dual edge lengths are positive (barycentric dual guarantee).
        const dual_lengths = mesh.dual_edge_volumes;
        for (dual_lengths) |dl| {
            try testing.expect(dl > 0.0);
        }

        // All dual vertex areas are positive.
        const dual_areas = mesh.vertices.slice().items(.dual_volume);
        for (dual_areas) |da| {
            try testing.expect(da > 0.0);
        }

        // Dual areas sum to total domain area.
        var total_dual: f64 = 0;
        for (dual_areas) |da| total_dual += da;
        try testing.expectApproxEqRel(width * height, total_dual, 1e-12);

        // Face areas sum to total domain area.
        var total_face: f64 = 0;
        for (areas) |a| total_face += a;
        try testing.expectApproxEqRel(width * height, total_face, 1e-12);
    }
}

test "cotangent Laplacian is robust on random grid dimensions (50 trials)" {
    // Stress test: Δ₀ on random meshes with random 0-forms must:
    //   1. Produce zero for constant functions (Δ₀(const) = 0)
    //   2. Be positive-semidefinite: ⟨ω, Δ₀ω⟩_★₀ ≥ 0
    const allocator = testing.allocator;
    const context_mod = @import("../operators/context.zig");
    const cochain_mod = @import("../forms/cochain.zig");

    var rng = std.Random.DefaultPrng.init(0xE5B_C07_00);

    for (0..50) |_| {
        const nx: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 15));
        const ny: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 15));
        const width = rng.random().float(f64) * 99.9 + 0.1;
        const height = rng.random().float(f64) * 99.9 + 0.1;

        var mesh = try Mesh(2, 2).uniform_grid(allocator, nx, ny, width, height);
        defer mesh.deinit(allocator);

        const operator_context = try context_mod.OperatorContext(Mesh(2, 2)).init(allocator, &mesh);
        defer operator_context.deinit();
        try operator_context.withLaplacian(0);

        const PrimalC0 = cochain_mod.Cochain(Mesh(2, 2), 0, cochain_mod.Primal);

        // Δ₀(constant) = 0.
        {
            var omega = try PrimalC0.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = 42.0;

            var result = try operator_context.laplacian(0).apply(allocator, omega);
            defer result.deinit(allocator);

            for (result.values) |v| {
                try testing.expectApproxEqAbs(@as(f64, 0.0), v, 1e-9);
            }
        }

        // Positive-semidefiniteness on a random 0-form.
        {
            var omega = try PrimalC0.init(allocator, &mesh);
            defer omega.deinit(allocator);
            for (omega.values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

            var lap_omega = try operator_context.laplacian(0).apply(allocator, omega);
            defer lap_omega.deinit(allocator);

            const dual_areas = mesh.vertices.slice().items(.dual_volume);
            var inner: f64 = 0;
            for (omega.values, lap_omega.values, dual_areas) |w, lw, area| {
                inner += w * lw * area;
            }
            try testing.expect(inner >= -1e-8);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Test helpers — minimal 3D mesh construction
// ═══════════════════════════════════════════════════════════════════════════

/// Build a single-tetrahedron Mesh(3, 3) for testing boundary operator properties.
///
/// Vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
/// Edges (6): all pairs of vertices, lexicographically ordered.
/// Faces (4): boundary faces of the tet, oriented so ∂₂∂₃ = 0 holds.
/// Tets (1): the single tetrahedron.
///
/// Orientation convention for a tet (v₀, v₁, v₂, v₃):
///   Face i (opposite vᵢ) gets sign (-1)^i in ∂₃.
///   Face vertices are ordered by skipping vᵢ from the sorted vertex list.
///   Each face's boundary edges in ∂₂ are oriented consistently with ∂₁.
fn build_single_tet(allocator: std.mem.Allocator) !Mesh(3, 3) {
    const M = Mesh(3, 3);

    // -- Vertices --
    var vertices = std.MultiArrayList(Vertex(M.embedding_dimension)){};
    try vertices.ensureTotalCapacity(allocator, 4);
    errdefer vertices.deinit(allocator);
    vertices.appendAssumeCapacity(.{ .coords = .{ 0, 0, 0 }, .dual_volume = 0.25 });
    vertices.appendAssumeCapacity(.{ .coords = .{ 1, 0, 0 }, .dual_volume = 0.25 });
    vertices.appendAssumeCapacity(.{ .coords = .{ 0, 1, 0 }, .dual_volume = 0.25 });
    vertices.appendAssumeCapacity(.{ .coords = .{ 0, 0, 1 }, .dual_volume = 0.25 });

    // -- Edges (6) --
    // Lexicographic: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    const edge_verts = [6][2]u32{
        .{ 0, 1 }, .{ 0, 2 }, .{ 0, 3 },
        .{ 1, 2 }, .{ 1, 3 }, .{ 2, 3 },
    };

    var edges = std.MultiArrayList(Simplex(M.embedding_dimension, M.topological_dimension, 1)){};
    try edges.ensureTotalCapacity(allocator, 6);
    errdefer edges.deinit(allocator);
    for (edge_verts) |ev| {
        const p0 = vertices.slice().items(.coords)[ev[0]];
        const p1 = vertices.slice().items(.coords)[ev[1]];
        edges.appendAssumeCapacity(.{
            .vertices = ev,
            .volume = M.euclidean_distance(p0, p1),
            .barycenter = M.point_midpoint(p0, p1),
        });
    }

    // Dual edge volumes (placeholder for 3D)
    const dual_edge_volumes = try allocator.alloc(f64, 6);
    errdefer allocator.free(dual_edge_volumes);
    @memset(dual_edge_volumes, 0.1);

    // -- ∂₁ (6 edges × 4 vertices) --
    // Each edge row: tail = -1, head = +1.
    var boundary_1: BoundaryMatrix = undefined;
    {
        var boundary_1_dense = try DenseBoundaryMatrix.init(allocator, 6, 4, 12);
        errdefer boundary_1_dense.deinit(allocator);
        for (0..6) |e| {
            boundary_1_dense.row_ptr[e] = @intCast(2 * e);
            boundary_1_dense.col_idx[2 * e] = edge_verts[e][0];
            boundary_1_dense.values[2 * e] = -1;
            boundary_1_dense.col_idx[2 * e + 1] = edge_verts[e][1];
            boundary_1_dense.values[2 * e + 1] = 1;
        }
        boundary_1_dense.row_ptr[6] = 12;
        boundary_1 = try BoundaryMatrix.fromBoundaryCsr(allocator, 1, boundary_1_dense);
        boundary_1_dense.deinit(allocator);
    }
    errdefer boundary_1.deinit(allocator);

    // -- Faces (4) --
    // Face i is opposite vertex i. Vertices of face i (sorted):
    //   face 0: (1,2,3), face 1: (0,2,3), face 2: (0,1,3), face 3: (0,1,2)
    const face_verts = [4][3]u32{
        .{ 1, 2, 3 }, // opposite v0
        .{ 0, 2, 3 }, // opposite v1
        .{ 0, 1, 3 }, // opposite v2
        .{ 0, 1, 2 }, // opposite v3
    };

    var faces = std.MultiArrayList(Simplex(M.embedding_dimension, M.topological_dimension, 2)){};
    try faces.ensureTotalCapacity(allocator, 4);
    errdefer faces.deinit(allocator);
    for (face_verts) |fv| {
        const p0 = vertices.slice().items(.coords)[fv[0]];
        const p1 = vertices.slice().items(.coords)[fv[1]];
        const p2 = vertices.slice().items(.coords)[fv[2]];
        faces.appendAssumeCapacity(.{
            .vertices = fv,
            .volume = M.triangle_area(p0, p1, p2),
            .barycenter = M.triangle_barycenter(p0, p1, p2),
        });
    }

    // -- ∂₂ (4 faces × 6 edges) --
    // For each face (a,b,c), the boundary edges are (a,b), (a,c), (b,c)
    // with signs from the orientation: +ab, -ac, +bc.
    //
    // Edge index lookup: edge (i,j) with i<j → index in lexicographic order.
    const edge_index = struct {
        fn f(i: u32, j: u32) u32 {
            // Lexicographic index for edge (i,j) among 4 vertices.
            // (0,1)=0, (0,2)=1, (0,3)=2, (1,2)=3, (1,3)=4, (2,3)=5
            std.debug.assert(i < j);
            return switch (i) {
                0 => j - 1,
                1 => j + 1,
                2 => 5,
                else => unreachable,
            };
        }
    }.f;

    var boundary_2: BoundaryMatrix = undefined;
    {
        var boundary_2_dense = try DenseBoundaryMatrix.init(allocator, 4, 6, 12);
        errdefer boundary_2_dense.deinit(allocator);
        for (0..4) |fi| {
            const fv = face_verts[fi];
            // Edges of face (a,b,c): (a,b) +1, (a,c) -1, (b,c) +1
            const e0 = edge_index(fv[0], fv[1]);
            const e1 = edge_index(fv[0], fv[2]);
            const e2 = edge_index(fv[1], fv[2]);

            // Store in sorted column order
            var cols: [3]u32 = .{ e0, e1, e2 };
            var vals: [3]i8 = .{ 1, -1, 1 };

            // Sort by column index (insertion sort on 3 elements)
            if (cols[0] > cols[1]) {
                std.mem.swap(u32, &cols[0], &cols[1]);
                std.mem.swap(i8, &vals[0], &vals[1]);
            }
            if (cols[1] > cols[2]) {
                std.mem.swap(u32, &cols[1], &cols[2]);
                std.mem.swap(i8, &vals[1], &vals[2]);
            }
            if (cols[0] > cols[1]) {
                std.mem.swap(u32, &cols[0], &cols[1]);
                std.mem.swap(i8, &vals[0], &vals[1]);
            }

            boundary_2_dense.row_ptr[fi] = @intCast(3 * fi);
            inline for (0..3) |j| {
                boundary_2_dense.col_idx[3 * fi + j] = cols[j];
                boundary_2_dense.values[3 * fi + j] = vals[j];
            }
        }
        boundary_2_dense.row_ptr[4] = 12;
        boundary_2 = try BoundaryMatrix.fromBoundaryCsr(allocator, 2, boundary_2_dense);
        boundary_2_dense.deinit(allocator);
    }
    errdefer boundary_2.deinit(allocator);

    // -- ∂₃ (1 tet × 4 faces) --
    // Sign for face i in ∂₃: (-1)^i.
    var boundary_3: BoundaryMatrix = undefined;
    {
        var boundary_3_dense = try DenseBoundaryMatrix.init(allocator, 1, 4, 4);
        errdefer boundary_3_dense.deinit(allocator);
        boundary_3_dense.row_ptr[0] = 0;
        boundary_3_dense.row_ptr[1] = 4;
        for (0..4) |i| {
            boundary_3_dense.col_idx[i] = @intCast(i);
            boundary_3_dense.values[i] = if (i % 2 == 0) @as(i8, 1) else @as(i8, -1);
        }
        boundary_3 = try BoundaryMatrix.fromBoundaryCsr(allocator, 3, boundary_3_dense);
        boundary_3_dense.deinit(allocator);
    }
    errdefer boundary_3.deinit(allocator);

    // -- Tets (1) --
    var tets = std.MultiArrayList(Simplex(M.embedding_dimension, M.topological_dimension, 3)){};
    try tets.ensureTotalCapacity(allocator, 1);
    errdefer tets.deinit(allocator);
    tets.appendAssumeCapacity(.{
        .vertices = .{ 0, 1, 2, 3 },
        .volume = 1.0 / 6.0,
        .barycenter = .{ 0.25, 0.25, 0.25 },
    });

    // Whitney mass and preconditioner are placeholders for 3D
    // (Whitney 1-form mass on tets is a different formula, tracked by #82).
    const dummy_mass = try sparse.CsrMatrix(f64).init(allocator, 0, 0, 0);
    errdefer dummy_mass.deinit(allocator);

    return .{
        .vertices = vertices,
        .simplex_lists = .{ edges, faces, tets },
        .boundaries = .{ boundary_1, boundary_2, boundary_3 },
        .dual_edge_volumes = dual_edge_volumes,
        .boundary_edges = &.{},
        .whitney_mass_1 = dummy_mass,
        .preconditioner_1 = &.{},
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// 3D mesh tests
// ═══════════════════════════════════════════════════════════════════════════

test "single tet entity counts" {
    const allocator = testing.allocator;
    var mesh = try build_single_tet(allocator);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(u32, 4), mesh.num_vertices());
    try testing.expectEqual(@as(u32, 6), mesh.num_edges());
    try testing.expectEqual(@as(u32, 4), mesh.num_faces());
    try testing.expectEqual(@as(u32, 1), mesh.num_tets());
}

test "uniform tetrahedral grid 1x1x1 entity counts" {
    const allocator = testing.allocator;
    var mesh = try Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(u32, 8), mesh.num_vertices());
    try testing.expectEqual(@as(u32, 19), mesh.num_edges());
    try testing.expectEqual(@as(u32, 18), mesh.num_faces());
    try testing.expectEqual(@as(u32, 6), mesh.num_tets());
}

test "uniform tetrahedral grid boundary of boundary is zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh(3, 3).uniform_tetrahedral_grid(allocator, 2, 2, 1, 2.0, 2.0, 1.0);
    defer mesh.deinit(allocator);

    var vertex_sum = try allocator.alloc(i32, mesh.num_vertices());
    defer allocator.free(vertex_sum);

    for (0..mesh.num_faces()) |f| {
        @memset(vertex_sum, 0);
        const face_row = mesh.boundary(2).row(@intCast(f));
        for (face_row.cols, 0..) |edge_idx, face_entry_idx| {
            const face_sign = face_row.sign(face_entry_idx);
            const edge_row = mesh.boundary(1).row(edge_idx);
            for (edge_row.cols, 0..) |vertex_idx, edge_entry_idx| {
                const edge_sign = edge_row.sign(edge_entry_idx);
                vertex_sum[vertex_idx] += @as(i32, face_sign) * @as(i32, edge_sign);
            }
        }
        for (vertex_sum) |sum| {
            try testing.expectEqual(@as(i32, 0), sum);
        }
    }

    var edge_sum = try allocator.alloc(i32, mesh.num_edges());
    defer allocator.free(edge_sum);

    for (0..mesh.num_tets()) |t| {
        @memset(edge_sum, 0);
        const tet_row = mesh.boundary(3).row(@intCast(t));
        for (tet_row.cols, 0..) |face_idx, tet_entry_idx| {
            const tet_sign = tet_row.sign(tet_entry_idx);
            const face_row = mesh.boundary(2).row(face_idx);
            for (face_row.cols, 0..) |edge_idx, face_entry_idx| {
                const face_sign = face_row.sign(face_entry_idx);
                edge_sum[edge_idx] += @as(i32, tet_sign) * @as(i32, face_sign);
            }
        }
        for (edge_sum) |sum| {
            try testing.expectEqual(@as(i32, 0), sum);
        }
    }
}

test "uniform tetrahedral grid Euler characteristic is one for random boxes" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x7E7_3D81);

    for (0..25) |_| {
        const nx: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 4));
        const ny: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 4));
        const nz: u32 = @intCast(rng.random().intRangeAtMost(u32, 1, 4));

        var mesh = try Mesh(3, 3).uniform_tetrahedral_grid(allocator, nx, ny, nz, 1.0, 1.0, 1.0);
        defer mesh.deinit(allocator);

        const chi = @as(i64, mesh.num_vertices()) -
            @as(i64, mesh.num_edges()) +
            @as(i64, mesh.num_faces()) -
            @as(i64, mesh.num_tets());
        try testing.expectEqual(@as(i64, 1), chi);
    }
}

test "uniform tetrahedral grid geometric measures are positive and conservative" {
    const allocator = testing.allocator;
    const width = 2.0;
    const height = 1.5;
    const depth = 0.75;

    var mesh = try Mesh(3, 3).uniform_tetrahedral_grid(allocator, 2, 3, 4, width, height, depth);
    defer mesh.deinit(allocator);

    const edge_lengths = mesh.simplices(1).items(.volume);
    for (edge_lengths) |length| {
        try testing.expect(length > 0.0);
    }

    const face_areas = mesh.simplices(2).items(.volume);
    for (face_areas) |area| {
        try testing.expect(area > 0.0);
    }

    const tet_volumes = mesh.simplices(3).items(.volume);
    var total_tet_volume: f64 = 0.0;
    for (tet_volumes) |tet_volume| {
        try testing.expect(tet_volume > 0.0);
        total_tet_volume += tet_volume;
    }
    try testing.expectApproxEqAbs(width * height * depth, total_tet_volume, 1e-12);

    const dual_vertex_volumes = mesh.vertices.slice().items(.dual_volume);
    var total_dual_vertex_volume: f64 = 0.0;
    for (dual_vertex_volumes) |dual_volume| {
        try testing.expect(dual_volume > 0.0);
        total_dual_vertex_volume += dual_volume;
    }
    try testing.expectApproxEqAbs(width * height * depth, total_dual_vertex_volume, 1e-12);

    for (mesh.dual_edge_volumes) |dual_face_area| {
        try testing.expect(dual_face_area > 0.0);
    }
}

test "∂₁∂₂ = 0 on single tetrahedron" {
    // For each face, applying ∂₁ to ∂₂(face) must yield zero at every vertex.
    // This is the boundary-of-boundary identity for edges → faces.
    const allocator = testing.allocator;
    var mesh = try build_single_tet(allocator);
    defer mesh.deinit(allocator);

    var vertex_sum = try allocator.alloc(i32, mesh.num_vertices());
    defer allocator.free(vertex_sum);

    for (0..mesh.num_faces()) |f| {
        @memset(vertex_sum, 0);
        const face_row = mesh.boundary(2).row(@intCast(f));
        for (face_row.cols, 0..) |edge_idx, face_entry_idx| {
            const face_sign = face_row.sign(face_entry_idx);
            const edge_row = mesh.boundary(1).row(edge_idx);
            for (edge_row.cols, 0..) |vert_idx, edge_entry_idx| {
                const edge_sign = edge_row.sign(edge_entry_idx);
                vertex_sum[vert_idx] += @as(i32, face_sign) * @as(i32, edge_sign);
            }
        }
        for (vertex_sum) |s| {
            try testing.expectEqual(@as(i32, 0), s);
        }
    }
}

test "∂₂∂₃ = 0 on single tetrahedron" {
    // For the single tet, applying ∂₂ to ∂₃(tet) must yield zero at every edge.
    // This is the boundary-of-boundary identity for faces → tets.
    const allocator = testing.allocator;
    var mesh = try build_single_tet(allocator);
    defer mesh.deinit(allocator);

    var edge_sum = try allocator.alloc(i32, mesh.num_edges());
    defer allocator.free(edge_sum);

    for (0..mesh.num_tets()) |t| {
        @memset(edge_sum, 0);
        const tet_row = mesh.boundary(3).row(@intCast(t));
        for (tet_row.cols, 0..) |face_idx, tet_entry_idx| {
            const tet_sign = tet_row.sign(tet_entry_idx);
            const face_row = mesh.boundary(2).row(face_idx);
            for (face_row.cols, 0..) |edge_idx, face_entry_idx| {
                const face_sign = face_row.sign(face_entry_idx);
                edge_sum[edge_idx] += @as(i32, tet_sign) * @as(i32, face_sign);
            }
        }
        for (edge_sum) |s| {
            try testing.expectEqual(@as(i32, 0), s);
        }
    }
}

test "uniform_grid 1×1 is the smallest valid grid" {
    // nx=0, ny=0, width=0, height≤0 all panic (precondition violations).
    // This test verifies the smallest valid grid constructs successfully.
    const allocator = testing.allocator;
    var mesh = try Mesh(2, 2).uniform_grid(allocator, 1, 1, 0.001, 0.001);
    defer mesh.deinit(allocator);
    try testing.expectEqual(@as(u32, 4), mesh.num_vertices());
    try testing.expectEqual(@as(u32, 2), mesh.num_faces());
}
