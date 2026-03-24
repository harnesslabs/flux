//! Simplicial mesh topology and geometry.
//!
//! Provides `Mesh(n)`, a 2D triangulation embedded in ℝⁿ with SoA entity
//! storage, oriented boundary operators ∂₁ and ∂₂ in CSR format, and
//! circumcentric dual geometry (dual edge lengths, dual vertex areas).

const std = @import("std");
const testing = std.testing;
const sparse = @import("../math/sparse.zig");

/// Boundary operators use i8-valued CSR matrices with entries in {−1, 0, +1}.
pub const BoundaryMatrix = sparse.CsrMatrix(i8);

// ───────────────────────────────────────────────────────────────────────────
// Mesh
// ───────────────────────────────────────────────────────────────────────────

/// Simplicial mesh parameterized on embedding dimension `n`.
///
/// Represents a 2D triangulation embedded in ℝⁿ. Topological entities
/// (vertices, edges, faces) use `std.MultiArrayList` for SoA cache-friendly
/// layout. Boundary operators ∂₁ and ∂₂ are stored in CSR format.
pub fn Mesh(comptime n: usize) type {
    comptime {
        if (n < 1) @compileError("Mesh dimension must be at least 1");
    }

    return struct {
        const Self = @This();

        /// Embedding dimension.
        pub const dimension = n;

        // -- Entity types stored in SoA layout via MultiArrayList --

        pub const Vertex = struct {
            coords: [n]f64,
            /// Area of the circumcentric dual cell (Voronoi region) around this vertex.
            dual_area: f64,
        };

        pub const Edge = struct {
            /// Vertex indices `[tail, head]` defining the edge orientation.
            vertices: [2]u32,
            /// Euclidean length of this edge.
            length: f64,
            /// Length of the circumcentric dual edge.
            /// For interior edges: distance between circumcenters of the two adjacent faces.
            /// For boundary edges: distance from the adjacent face's circumcenter to the edge midpoint.
            dual_length: f64,
        };

        /// A triangular face with CCW-oriented vertex indices and geometric data.
        pub const Face = struct {
            /// Vertex indices `[v0, v1, v2]` in counter-clockwise orientation.
            vertices: [3]u32,
            /// Area of this triangle.
            area: f64,
            /// Circumcenter coordinates.
            circumcenter: [n]f64,
        };

        // -- Storage --

        vertices: std.MultiArrayList(Vertex),
        edges: std.MultiArrayList(Edge),
        faces: std.MultiArrayList(Face),

        /// ∂₁: `n_edges × n_vertices`. Row `e` has nonzeros at tail (−1) and head (+1).
        boundary_1: BoundaryMatrix,
        /// ∂₂: `n_faces × n_edges`. Row `f` has 3 nonzeros for the oriented boundary edges.
        boundary_2: BoundaryMatrix,

        /// Indices of edges on the mesh boundary (adjacent to exactly one face).
        /// Precomputed during construction for use by boundary condition routines.
        boundary_edges: []u32,

        // -- Lifetime --

        /// Free all entity storage and boundary matrices.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.boundary_edges);
            self.vertices.deinit(allocator);
            self.edges.deinit(allocator);
            self.faces.deinit(allocator);
            self.boundary_1.deinit(allocator);
            self.boundary_2.deinit(allocator);
        }

        // -- Accessors --

        /// The topological dimension of this mesh (always 2 for now).
        pub const topological_dimension = 2;

        /// Return the boundary operator ∂ₖ for the given degree.
        ///
        /// ∂₁ maps edges → vertices, ∂₂ maps faces → edges. Stored in
        /// coboundary orientation (rows indexed by higher-dimensional cells),
        /// so ∂ₖ also serves as the exterior derivative dₖ₋₁.
        pub fn boundary(self: Self, comptime k: comptime_int) BoundaryMatrix {
            return switch (k) {
                1 => self.boundary_1,
                2 => self.boundary_2,
                else => @compileError(std.fmt.comptimePrint(
                    "no boundary operator ∂_{d} on a {d}-dimensional mesh",
                    .{ k, topological_dimension },
                )),
            };
        }

        /// Number of vertices in the mesh.
        pub fn num_vertices(self: Self) u32 {
            return @intCast(self.vertices.len);
        }

        /// Number of edges in the mesh.
        pub fn num_edges(self: Self) u32 {
            return @intCast(self.edges.len);
        }

        /// Number of faces in the mesh.
        pub fn num_faces(self: Self) u32 {
            return @intCast(self.faces.len);
        }

        // ───────────────────────────────────────────────────────────────────
        // Uniform grid constructor
        // ───────────────────────────────────────────────────────────────────

        /// Construct a uniform triangulated rectangular grid on `[0, width] × [0, height]`.
        ///
        /// The rectangle is divided into `nx × ny` rectangular cells, each split
        /// by a SW-to-NE diagonal into two CCW-oriented triangles (`2 * nx * ny` total).
        /// All geometric data — lengths, areas, circumcenters, dual volumes — is
        /// computed during construction.
        ///
        /// **Note on circumcentric dual degeneracy:** both triangles in each cell
        /// share the same circumcenter (center of the cell) because they are right
        /// triangles whose hypotenuse is the cell diagonal. The diagonal edges
        /// therefore have dual length 0. This is documented and expected; the Hodge
        /// star in M2 will need to handle this boundary case.
        pub fn uniform_grid(
            allocator: std.mem.Allocator,
            nx: u32,
            ny: u32,
            width: f64,
            height: f64,
        ) !Self {
            std.debug.assert(nx > 0);
            std.debug.assert(ny > 0);
            std.debug.assert(width > 0);
            std.debug.assert(height > 0);

            const dx = width / @as(f64, @floatFromInt(nx));
            const dy = height / @as(f64, @floatFromInt(ny));

            const vertex_count: u32 = (nx + 1) * (ny + 1);
            const horizontal_edge_count: u32 = nx * (ny + 1);
            const vertical_edge_count: u32 = (nx + 1) * ny;
            const diagonal_edge_count: u32 = nx * ny;
            const edge_count: u32 = horizontal_edge_count + vertical_edge_count + diagonal_edge_count;
            const face_count: u32 = 2 * nx * ny;

            // -- Allocate entity storage --

            var vertices = std.MultiArrayList(Vertex){};
            try vertices.ensureTotalCapacity(allocator, vertex_count);
            errdefer vertices.deinit(allocator);

            var edges_list = std.MultiArrayList(Edge){};
            try edges_list.ensureTotalCapacity(allocator, edge_count);
            errdefer edges_list.deinit(allocator);

            var faces_list = std.MultiArrayList(Face){};
            try faces_list.ensureTotalCapacity(allocator, face_count);
            errdefer faces_list.deinit(allocator);

            // -- Populate vertices --
            // Indexed as vertex_index(i, j) = i * (ny + 1) + j.
            for (0..nx + 1) |i_u| {
                const fi: f64 = @floatFromInt(i_u);
                for (0..ny + 1) |j_u| {
                    const fj: f64 = @floatFromInt(j_u);
                    var coords: [n]f64 = @splat(0);
                    coords[0] = fi * dx;
                    coords[1] = fj * dy;
                    vertices.appendAssumeCapacity(.{ .coords = coords, .dual_area = 0 });
                }
            }

            // -- Populate edges --
            // Horizontal: h(i,j) = j*nx + i, tail = v(i,j), head = v(i+1,j)
            for (0..ny + 1) |j_u| {
                for (0..nx) |i_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i + 1, j, ny) },
                        .length = 0,
                        .dual_length = 0,
                    });
                }
            }
            // Vertical: vert(i,j) = horizontal_edge_count + i*ny + j
            for (0..nx + 1) |i_u| {
                for (0..ny) |j_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i, j + 1, ny) },
                        .length = 0,
                        .dual_length = 0,
                    });
                }
            }
            // Diagonal (SW→NE): d(i,j) = horizontal_edge_count + vertical_edge_count + i*ny + j
            for (0..nx) |i_u| {
                for (0..ny) |j_u| {
                    const i: u32 = @intCast(i_u);
                    const j: u32 = @intCast(j_u);
                    edges_list.appendAssumeCapacity(.{
                        .vertices = .{ vertex_index(i, j, ny), vertex_index(i + 1, j + 1, ny) },
                        .length = 0,
                        .dual_length = 0,
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

                    var zero_cc: [n]f64 = @splat(0);
                    _ = &zero_cc;

                    // Lower-right triangle: SW → SE → NE
                    faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, se, ne }, .area = 0, .circumcenter = zero_cc });
                    // Upper-left triangle:  SW → NE → NW
                    faces_list.appendAssumeCapacity(.{ .vertices = .{ sw, ne, nw }, .area = 0, .circumcenter = zero_cc });
                }
            }

            // -- Build ∂₁ (n_edges × n_vertices) --
            // Each edge row has exactly 2 nonzeros: tail = −1, head = +1.
            var boundary_1 = try BoundaryMatrix.init(allocator, edge_count, vertex_count, 2 * edge_count);
            errdefer boundary_1.deinit(allocator);
            {
                const edge_verts = edges_list.slice().items(.vertices);
                for (0..edge_count) |e| {
                    boundary_1.row_ptr[e] = @intCast(2 * e);
                    // tail < head for the uniform grid, so columns are sorted
                    boundary_1.col_idx[2 * e] = edge_verts[e][0];
                    boundary_1.values[2 * e] = -1;
                    boundary_1.col_idx[2 * e + 1] = edge_verts[e][1];
                    boundary_1.values[2 * e + 1] = 1;
                }
                boundary_1.row_ptr[edge_count] = 2 * edge_count;
            }

            // -- Build ∂₂ (n_faces × n_edges) --
            // Each face row has exactly 3 nonzeros for the oriented boundary edges.
            //
            // Lower-right (SW,SE,NE): +h(i,j), +vert(i+1,j), −diag(i,j)
            // Upper-left  (SW,NE,NW): −h(i,j+1), −vert(i,j), +diag(i,j)
            //
            // Column ordering: horizontal < vertical < diagonal (always sorted).
            var boundary_2 = try BoundaryMatrix.init(allocator, face_count, edge_count, 3 * face_count);
            errdefer boundary_2.deinit(allocator);
            {
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
                        boundary_2.row_ptr[f_idx] = 3 * f_idx;
                        boundary_2.col_idx[3 * f_idx + 0] = h_ij;
                        boundary_2.values[3 * f_idx + 0] = 1;
                        boundary_2.col_idx[3 * f_idx + 1] = v_ip1_j;
                        boundary_2.values[3 * f_idx + 1] = 1;
                        boundary_2.col_idx[3 * f_idx + 2] = d_ij;
                        boundary_2.values[3 * f_idx + 2] = -1;
                        f_idx += 1;

                        // Upper-left triangle
                        boundary_2.row_ptr[f_idx] = 3 * f_idx;
                        boundary_2.col_idx[3 * f_idx + 0] = h_i_jp1;
                        boundary_2.values[3 * f_idx + 0] = -1;
                        boundary_2.col_idx[3 * f_idx + 1] = v_i_j;
                        boundary_2.values[3 * f_idx + 1] = -1;
                        boundary_2.col_idx[3 * f_idx + 2] = d_ij;
                        boundary_2.values[3 * f_idx + 2] = 1;
                        f_idx += 1;
                    }
                }
                boundary_2.row_ptr[face_count] = 3 * face_count;
            }

            // -- Compute primal geometry --

            const coords = vertices.slice().items(.coords);

            // Edge lengths
            {
                const edge_verts = edges_list.slice().items(.vertices);
                const lengths = edges_list.slice().items(.length);
                for (0..edge_count) |e| {
                    lengths[e] = euclidean_distance(coords[edge_verts[e][0]], coords[edge_verts[e][1]]);
                }
            }

            // Face areas and circumcenters
            {
                const face_verts = faces_list.slice().items(.vertices);
                const areas = faces_list.slice().items(.area);
                const circumcenters = faces_list.slice().items(.circumcenter);
                for (0..face_count) |f| {
                    const p0 = coords[face_verts[f][0]];
                    const p1 = coords[face_verts[f][1]];
                    const p2 = coords[face_verts[f][2]];
                    areas[f] = triangle_area(p0, p1, p2);
                    circumcenters[f] = triangle_circumcenter(p0, p1, p2);
                }
            }

            // -- Compute circumcentric dual geometry --

            // Dual edge lengths require edge→face adjacency.
            // Each edge borders at most 2 faces; boundary edges border exactly 1.
            var boundary_edge_buf: []u32 = &.{};
            {
                var edge_face_count = try allocator.alloc(u8, edge_count);
                defer allocator.free(edge_face_count);
                var edge_face_0 = try allocator.alloc(u32, edge_count);
                defer allocator.free(edge_face_0);
                var edge_face_1 = try allocator.alloc(u32, edge_count);
                defer allocator.free(edge_face_1);
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

                const circumcenters = faces_list.slice().items(.circumcenter);
                const edge_verts = edges_list.slice().items(.vertices);
                const dual_lengths = edges_list.slice().items(.dual_length);

                // Count boundary edges (adjacent to exactly one face) and
                // compute dual lengths in one pass.
                var boundary_count: u32 = 0;
                for (0..edge_count) |e| {
                    if (edge_face_count[e] == 2) {
                        dual_lengths[e] = euclidean_distance(
                            circumcenters[edge_face_0[e]],
                            circumcenters[edge_face_1[e]],
                        );
                    } else if (edge_face_count[e] == 1) {
                        const mid = point_midpoint(coords[edge_verts[e][0]], coords[edge_verts[e][1]]);
                        dual_lengths[e] = euclidean_distance(circumcenters[edge_face_0[e]], mid);
                        boundary_count += 1;
                    }
                    // edge_face_count == 0 should not occur in a valid mesh; dual_length stays 0
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

            // Dual vertex areas via the cotangent formula.
            //
            // For triangle (v₀, v₁, v₂) with edge lengths l₀₁, l₁₂, l₂₀:
            //   cot(αᵢ) = (lⱼₖ² + lₖᵢ² − lᵢⱼ²) / (4 · area)
            //
            //   dual_area[vᵢ] += (cot(αⱼ) · lᵢₖ² + cot(αₖ) · lᵢⱼ²) / 8
            //
            // This gives the circumcentric (Voronoi) dual area contribution.
            // Reference: Meyer et al., "Discrete Differential-Geometry Operators
            // for Triangulated 2-Manifolds" (2002), §4.
            {
                const dual_areas = vertices.slice().items(.dual_area);
                @memset(dual_areas, 0);

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
                    std.debug.assert(area_4 > 0); // degenerate triangle

                    const cot0 = (l01_sq + l20_sq - l12_sq) / area_4;
                    const cot1 = (l01_sq + l12_sq - l20_sq) / area_4;
                    const cot2 = (l12_sq + l20_sq - l01_sq) / area_4;

                    dual_areas[vs[0]] += (cot1 * l20_sq + cot2 * l01_sq) / 8.0;
                    dual_areas[vs[1]] += (cot2 * l01_sq + cot0 * l12_sq) / 8.0;
                    dual_areas[vs[2]] += (cot0 * l12_sq + cot1 * l20_sq) / 8.0;
                }
            }

            return Self{
                .vertices = vertices,
                .edges = edges_list,
                .faces = faces_list,
                .boundary_1 = boundary_1,
                .boundary_2 = boundary_2,
                .boundary_edges = boundary_edge_buf,
            };
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

        fn distance_squared(a: [n]f64, b: [n]f64) f64 {
            var sum: f64 = 0;
            inline for (0..n) |d| {
                const diff = b[d] - a[d];
                sum += diff * diff;
            }
            return sum;
        }

        fn euclidean_distance(a: [n]f64, b: [n]f64) f64 {
            return @sqrt(distance_squared(a, b));
        }

        fn point_midpoint(a: [n]f64, b: [n]f64) [n]f64 {
            var result: [n]f64 = undefined;
            inline for (0..n) |d| {
                result[d] = 0.5 * (a[d] + b[d]);
            }
            return result;
        }

        /// Signed area of a 2D triangle (positive for CCW orientation).
        /// Uses only the first two coordinate components.
        fn signed_triangle_area(p0: [n]f64, p1: [n]f64, p2: [n]f64) f64 {
            const dx1 = p1[0] - p0[0];
            const dy1 = p1[1] - p0[1];
            const dx2 = p2[0] - p0[0];
            const dy2 = p2[1] - p0[1];
            return (dx1 * dy2 - dx2 * dy1) / 2.0;
        }

        /// Unsigned area of a 2D triangle. Uses only the first two coordinate components.
        fn triangle_area(p0: [n]f64, p1: [n]f64, p2: [n]f64) f64 {
            return @abs(signed_triangle_area(p0, p1, p2));
        }

        /// Circumcenter of a 2D triangle. Uses only the first two coordinate components.
        fn triangle_circumcenter(a: [n]f64, b: [n]f64, c: [n]f64) [n]f64 {
            const ax = a[0];
            const ay = a[1];
            const bx = b[0];
            const by = b[1];
            const cx = c[0];
            const cy = c[1];

            const d_val = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
            std.debug.assert(@abs(d_val) > 1e-15); // degenerate triangle guard

            const a_sq = ax * ax + ay * ay;
            const b_sq = bx * bx + by * by;
            const c_sq = cx * cx + cy * cy;

            var result: [n]f64 = @splat(0);
            result[0] = (a_sq * (by - cy) + b_sq * (cy - ay) + c_sq * (ay - by)) / d_val;
            result[1] = (a_sq * (cx - bx) + b_sq * (ax - cx) + c_sq * (bx - ax)) / d_val;
            return result;
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

test "uniform grid entity counts" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2).uniform_grid(allocator, 3, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(u32, 20), mesh.num_vertices()); // (3+1)*(4+1)
    try testing.expectEqual(@as(u32, 43), mesh.num_edges()); // 3*5 + 4*4 + 3*4
    try testing.expectEqual(@as(u32, 24), mesh.num_faces()); // 2*3*4
}

test "boundary operator ∂₁ has exactly 2 nonzeros per row" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    for (0..mesh.num_edges()) |e| {
        const r = mesh.boundary_1.row(@intCast(e));
        try testing.expectEqual(@as(usize, 2), r.cols.len);
        // tail = −1, head = +1
        try testing.expectEqual(@as(i8, -1), r.vals[0]);
        try testing.expectEqual(@as(i8, 1), r.vals[1]);
        // columns sorted
        try testing.expect(r.cols[0] < r.cols[1]);
    }
}

test "boundary operator ∂₂ has exactly 3 nonzeros per row" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    for (0..mesh.num_faces()) |f| {
        const r = mesh.boundary_2.row(@intCast(f));
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
        var mesh = try Mesh(2).uniform_grid(allocator, size[0], size[1], 1.0, 1.0);
        defer mesh.deinit(allocator);

        var vertex_sum = try allocator.alloc(i32, mesh.num_vertices());
        defer allocator.free(vertex_sum);

        for (0..mesh.num_faces()) |f| {
            @memset(vertex_sum, 0);

            const face_row = mesh.boundary_2.row(@intCast(f));
            for (face_row.cols, face_row.vals) |edge_idx, face_sign| {
                const edge_row = mesh.boundary_1.row(edge_idx);
                for (edge_row.cols, edge_row.vals) |vert_idx, edge_sign| {
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
    var mesh = try Mesh(2).uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const lengths = mesh.edges.slice().items(.length);

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
    var mesh = try Mesh(2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const areas = mesh.faces.slice().items(.area);
    const expected_area: f64 = 0.5 * 0.5 / 2.0; // dx * dy / 2

    for (0..mesh.num_faces()) |f| {
        try testing.expectApproxEqAbs(expected_area, areas[f], 1e-15);
    }
}

test "circumcenters of right triangles are at hypotenuse midpoints" {
    const allocator = testing.allocator;
    var mesh = try Mesh(2).uniform_grid(allocator, 1, 1, 2.0, 2.0);
    defer mesh.deinit(allocator);

    const circumcenters = mesh.faces.slice().items(.circumcenter);

    // Single cell: both triangles have the same circumcenter at (1, 1),
    // the midpoint of the shared hypotenuse (which is also the cell center).
    for (0..mesh.num_faces()) |f| {
        try testing.expectApproxEqAbs(@as(f64, 1.0), circumcenters[f][0], 1e-15);
        try testing.expectApproxEqAbs(@as(f64, 1.0), circumcenters[f][1], 1e-15);
    }
}

test "dual vertex areas sum to total mesh area" {
    const allocator = testing.allocator;
    const width: f64 = 3.0;
    const height: f64 = 2.0;

    var mesh = try Mesh(2).uniform_grid(allocator, 5, 4, width, height);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_area);
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

    var mesh = try Mesh(2).uniform_grid(allocator, nx, ny, width, height);
    defer mesh.deinit(allocator);

    const dual_areas = mesh.vertices.slice().items(.dual_area);

    // Interior vertices should each have dual area = dx * dy
    for (1..nx) |i| {
        for (1..ny) |j| {
            const idx = Mesh(2).vertex_index(@intCast(i), @intCast(j), ny);
            try testing.expectApproxEqAbs(dx * dy, dual_areas[idx], 1e-13);
        }
    }
}

test "diagonal edges have zero dual length" {
    // Both triangles sharing a diagonal have the same circumcenter,
    // so the dual edge connecting those circumcenters has length 0.
    const allocator = testing.allocator;
    var mesh = try Mesh(2).uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const dual_lengths = mesh.edges.slice().items(.dual_length);

    // Diagonal edges start at index: horizontal_count + vertical_count
    const n_h: u32 = 3 * 4; // nx * (ny + 1) = 3 * 4 = 12
    const n_v: u32 = 4 * 3; // (nx + 1) * ny = 4 * 3 = 12
    const diag_start = n_h + n_v;

    for (diag_start..mesh.num_edges()) |e| {
        try testing.expectApproxEqAbs(@as(f64, 0), dual_lengths[e], 1e-15);
    }
}

test "Mesh(2) compiles at dimension 2" {
    // Compile-time check: Mesh(2) is a valid type.
    const M = Mesh(2);
    try testing.expect(M.dimension == 2);
}

test "Mesh(3) compiles at dimension 3" {
    // Compile-time check: Mesh(3) is a valid type (future-proofing).
    const M = Mesh(3);
    try testing.expect(M.dimension == 3);
}
