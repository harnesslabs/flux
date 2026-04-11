//! Benchmark suite for flux operators.
//!
//! Measures wall-clock performance of key operators on a moderately-sized mesh.
//! Reports median timings (robust to outliers) in both human-readable table
//! and machine-readable JSON formats.
//!
//! Usage:
//!   zig build bench                    — run default benchmarks, print results
//!   zig build bench -- --check         — compare default suite against baselines, exit 1 on regression
//!   zig build bench -- --update        — run default suite and overwrite baselines.json
//!   zig build bench -- --json          — print JSON to stdout (default: table to stderr)
//!   zig build bench -- --compare       — run same-run implementation comparisons

const std = @import("std");
const builtin = @import("builtin");
const flux = @import("flux");
const maxwell = @import("maxwell_example");
const diffusion_sphere = @import("diffusion_sphere_example");
/// Version of the benchmark result file schema.
const benchmark_suite_version: u32 = 4;
/// Version for operator and end-to-end rows after switching to calibrated samples.
const stable_benchmark_method_version: u32 = 2;
/// Version for tiny cochain microbenchmarks after switching to calibrated samples.
const micro_benchmark_method_version: u32 = 3;
/// Version for arithmetic scalar-vs-default comparisons after switching to interleaved same-run samples.
const arithmetic_benchmark_method_version: u32 = 3;
/// Version for incidence compare rows after switching to interleaved same-run samples.
const incidence_compare_method_version: u32 = 2;

const Mesh2D = flux.topology.Mesh(2, 2);
const SurfaceMesh = flux.topology.Mesh(3, 2);
const PrimalC0 = flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal);
const PrimalC1 = flux.forms.Cochain(Mesh2D, 1, flux.forms.Primal);
const PrimalC2 = flux.forms.Cochain(Mesh2D, 2, flux.forms.Primal);
const OperatorContext2D = flux.operators.context.OperatorContext(Mesh2D);
const SurfaceOperatorContext = flux.operators.context.OperatorContext(SurfaceMesh);
const MaxwellState2D = maxwell.State(2);
const ArithmeticMesh = struct {
    pub const embedding_dimension = 2;
    pub const topological_dimension = 2;

    vertices: u32,
    edges: u32,
    faces: u32,

    pub fn num_vertices(self: *const @This()) u32 {
        return self.vertices;
    }

    pub fn num_edges(self: *const @This()) u32 {
        return self.edges;
    }

    pub fn num_faces(self: *const @This()) u32 {
        return self.faces;
    }

    pub fn num_cells(self: *const @This(), comptime k: comptime_int) u32 {
        return switch (k) {
            0 => self.num_vertices(),
            1 => self.num_edges(),
            2 => self.num_faces(),
            else => @compileError("ArithmeticMesh only supports cell counts for degrees 0, 1, and 2"),
        };
    }
};
const ArithmeticC1 = flux.forms.Cochain(ArithmeticMesh, 1, flux.forms.Primal);

// ── Configuration ───────────────────────────────────────────────────────

const grid_nx: u32 = 50;
const grid_ny: u32 = 50;
const grid_width: f64 = 50.0;
const grid_height: f64 = 50.0;
const warmup_iterations: u32 = 10;
const measured_iterations: u32 = 100;
const arithmetic_scale_factor: f64 = 1.000001;
const arithmetic_case_lengths = [_]u32{ 1024, 16 * 1024, 128 * 1024 };
const arithmetic_case_labels = [_][]const u8{ "1k", "16k", "128k" };
const cavity_grid: u32 = 256;
const cavity_domain: f64 = 1.0;
const cavity_courant: f64 = 0.1;
const surface_patch_nx: u32 = 64;
const surface_patch_ny: u32 = 64;
const diffusion_surface_refinement: u32 = 5;
const diffusion_surface_final_time: f64 = 0.05;
const diffusion_surface_steps: u32 = 1;
const small_cochain_repetitions: u32 = 64;
const incidence_grid_nx: u32 = 1000;
const incidence_grid_ny: u32 = 1000;
const incidence_warmup_iterations: u32 = 5;
const incidence_measured_iterations: u32 = 20;
const target_sample_duration_ns: u64 = 5 * std.time.ns_per_ms;
const max_sample_repetitions: u32 = 1 << 20;
/// Maximum allowed regression before CI fails (0.20 = 20%).
const regression_threshold: f64 = 0.20;

// ── Benchmark definition ────────────────────────────────────────────────

const BenchmarkFn = *const fn (*BenchmarkContext) void;

const BenchmarkClass = enum {
    gate,
    info,
};

const BenchmarkScope = enum {
    public_library,
    shipped_example,
    compare_experiment,
};

const BenchmarkMode = enum {
    default,
    compare,
};

const BenchmarkDef = struct {
    name: []const u8,
    run: BenchmarkFn,
    minimum_repetitions: u32 = 1,
    version: u32 = stable_benchmark_method_version,
    class: BenchmarkClass = .gate,
    scope: BenchmarkScope,
};

const BenchmarkResult = struct {
    name: []const u8,
    version: u32,
    class: BenchmarkClass,
    repetitions: u32,
    median_ns: u64,
    min_ns: u64,
    max_ns: u64,
    iterations: u32,
};

const ArithmeticOperation = enum {
    add,
    scale,
    negate,
    inner_product,
};

const ArithmeticCompareDef = struct {
    name: []const u8,
    operation: ArithmeticOperation,
    case_index: usize,
    version: u32 = arithmetic_benchmark_method_version,
    scope: BenchmarkScope = .compare_experiment,
};

const ArithmeticComparisonResult = struct {
    name: []const u8,
    operation: ArithmeticOperation,
    case_label: []const u8,
    version: u32,
    scalar_median_ns: u64,
    default_median_ns: u64,
};

const IncidenceComparisonResult = struct {
    name: []const u8,
    version: u32 = incidence_compare_method_version,
    nnz: u32,
    dense_value_bytes: usize,
    generic_packed_value_bytes: usize,
    specialized_value_bytes: usize,
    dense_median_ns: u64,
    generic_packed_median_ns: u64,
    specialized_median_ns: u64,
};

const BenchmarkArgs = struct {
    mode: BenchmarkMode = .default,
    check_mode: bool = false,
    update_mode: bool = false,
    json_mode: bool = false,
};

const ArithmeticCase = struct {
    mesh: ArithmeticMesh,
    lhs: ArithmeticC1,
    rhs: ArithmeticC1,

    fn init(allocator: std.mem.Allocator, len: u32) !ArithmeticCase {
        var mesh = ArithmeticMesh{
            .vertices = len,
            .edges = len,
            .faces = len,
        };
        var lhs = try ArithmeticC1.init(allocator, &mesh);
        errdefer lhs.deinit(allocator);
        var rhs = try ArithmeticC1.init(allocator, &mesh);
        errdefer rhs.deinit(allocator);

        var rng = std.Random.DefaultPrng.init(@as(u64, len) * 0x9E37_79B9);
        fillRandom(&rng, lhs.values);
        fillRandom(&rng, rhs.values);

        return .{
            .mesh = mesh,
            .lhs = lhs,
            .rhs = rhs,
        };
    }

    fn deinit(self: *ArithmeticCase, allocator: std.mem.Allocator) void {
        self.rhs.deinit(allocator);
        self.lhs.deinit(allocator);
    }
};

const BenchmarkContext = struct {
    mesh: *Mesh2D,
    allocator: std.mem.Allocator,
    // Pre-allocated cochains for benchmarks that need them.
    // Benchmarks should not measure allocation time — only operator application.
    c0: PrimalC0,
    c1: PrimalC1,
    c2: PrimalC2,
    // Secondary cochains for binary operations.
    c0_other: PrimalC0,
    c1_other: PrimalC1,
    arithmetic_cases: [arithmetic_case_lengths.len]ArithmeticCase,
    cavity_mesh: *Mesh2D,
    cavity_state: MaxwellState2D,
    operator_context: *OperatorContext2D,
    surface_mesh: *SurfaceMesh,

    fn init(allocator: std.mem.Allocator, mesh: *Mesh2D) !BenchmarkContext {
        var c0 = try PrimalC0.init(allocator, mesh);
        errdefer c0.deinit(allocator);
        var c1 = try PrimalC1.init(allocator, mesh);
        errdefer c1.deinit(allocator);
        var c2 = try PrimalC2.init(allocator, mesh);
        errdefer c2.deinit(allocator);
        var c0_other = try PrimalC0.init(allocator, mesh);
        errdefer c0_other.deinit(allocator);
        var c1_other = try PrimalC1.init(allocator, mesh);
        errdefer c1_other.deinit(allocator);
        var arithmetic_cases: [arithmetic_case_lengths.len]ArithmeticCase = undefined;
        var arithmetic_cases_initialized: usize = 0;
        errdefer for (arithmetic_cases[0..arithmetic_cases_initialized]) |*arithmetic_case| {
            arithmetic_case.deinit(allocator);
        };
        inline for (arithmetic_case_lengths, 0..) |len, i| {
            arithmetic_cases[i] = try ArithmeticCase.init(allocator, len);
            arithmetic_cases_initialized = i + 1;
        }
        const cavity_mesh = try allocator.create(Mesh2D);
        errdefer allocator.destroy(cavity_mesh);
        cavity_mesh.* = try Mesh2D.uniform_grid(allocator, cavity_grid, cavity_grid, cavity_domain, cavity_domain);
        errdefer cavity_mesh.deinit(allocator);
        var cavity_state = try MaxwellState2D.init(allocator, cavity_mesh);
        errdefer cavity_state.deinit(allocator);
        try maxwell.seedReferenceMode(2, allocator, &cavity_state, cavityDt(), cavity_domain, cavity_domain);
        const operator_context = try OperatorContext2D.init(allocator, mesh);
        errdefer operator_context.deinit();
        const surface_mesh = try allocator.create(SurfaceMesh);
        errdefer allocator.destroy(surface_mesh);
        surface_mesh.* = try buildEmbeddedSurfacePatch(allocator, surface_patch_nx, surface_patch_ny);
        errdefer surface_mesh.deinit(allocator);
        _ = try operator_context.exteriorDerivative(flux.forms.Primal, 0);
        _ = try operator_context.exteriorDerivative(flux.forms.Primal, 1);
        _ = try operator_context.hodgeStar(0);
        _ = try operator_context.hodgeStar(1);
        _ = try operator_context.hodgeStar(2);
        _ = try operator_context.hodgeStarInverse(1);
        _ = try operator_context.laplacian(0);

        // Fill with deterministic pseudo-random values.
        var rng = std.Random.DefaultPrng.init(0xBEEF_CAFE);
        fillRandom(&rng, c0.values);
        fillRandom(&rng, c1.values);
        fillRandom(&rng, c2.values);
        fillRandom(&rng, c0_other.values);
        fillRandom(&rng, c1_other.values);

        return .{
            .mesh = mesh,
            .allocator = allocator,
            .c0 = c0,
            .c1 = c1,
            .c2 = c2,
            .c0_other = c0_other,
            .c1_other = c1_other,
            .arithmetic_cases = arithmetic_cases,
            .cavity_mesh = cavity_mesh,
            .cavity_state = cavity_state,
            .operator_context = operator_context,
            .surface_mesh = surface_mesh,
        };
    }

    fn deinit(self: *BenchmarkContext) void {
        self.surface_mesh.deinit(self.allocator);
        self.allocator.destroy(self.surface_mesh);
        self.operator_context.deinit();
        self.cavity_state.deinit(self.allocator);
        self.cavity_mesh.deinit(self.allocator);
        self.allocator.destroy(self.cavity_mesh);
        for (&self.arithmetic_cases) |*arithmetic_case| {
            arithmetic_case.deinit(self.allocator);
        }
        self.c1_other.deinit(self.allocator);
        self.c0_other.deinit(self.allocator);
        self.c2.deinit(self.allocator);
        self.c1.deinit(self.allocator);
        self.c0.deinit(self.allocator);
    }
};

fn fillRandom(rng: *std.Random.DefaultPrng, values: []f64) void {
    for (values) |*v| {
        v.* = rng.random().float(f64) * 200.0 - 100.0;
    }
}

fn cavityDt() f64 {
    return cavity_courant * (cavity_domain / @as(f64, @floatFromInt(cavity_grid)));
}

fn buildEmbeddedSurfacePatch(allocator: std.mem.Allocator, nx: u32, ny: u32) !SurfaceMesh {
    std.debug.assert(nx > 0);
    std.debug.assert(ny > 0);

    const vertex_count = (nx + 1) * (ny + 1);
    const face_count = 2 * nx * ny;
    const inv_sqrt2 = 0.7071067811865475;

    const vertices = try allocator.alloc([3]f64, vertex_count);
    defer allocator.free(vertices);
    const faces = try allocator.alloc([3]u32, face_count);
    defer allocator.free(faces);

    var vertex_write: usize = 0;
    for (0..nx + 1) |i_u| {
        const x = @as(f64, @floatFromInt(i_u)) / @as(f64, @floatFromInt(nx));
        for (0..ny + 1) |j_u| {
            const y = @as(f64, @floatFromInt(j_u)) / @as(f64, @floatFromInt(ny));
            vertices[vertex_write] = .{
                (x + y) * inv_sqrt2,
                (x - y) * inv_sqrt2,
                0.5 * y,
            };
            vertex_write += 1;
        }
    }
    std.debug.assert(vertex_write == vertex_count);

    var face_write: usize = 0;
    for (0..nx) |i_u| {
        const i: u32 = @intCast(i_u);
        for (0..ny) |j_u| {
            const j: u32 = @intCast(j_u);
            const sw = i * (ny + 1) + j;
            const se = (i + 1) * (ny + 1) + j;
            const nw = i * (ny + 1) + (j + 1);
            const ne = (i + 1) * (ny + 1) + (j + 1);
            faces[face_write] = .{ sw, se, ne };
            faces[face_write + 1] = .{ sw, ne, nw };
            face_write += 2;
        }
    }
    std.debug.assert(face_write == face_count);

    return SurfaceMesh.from_triangles(allocator, vertices, faces);
}

fn addScalarLoop(lhs: []f64, rhs: []const f64) void {
    for (lhs, rhs) |*left, right| {
        left.* += right;
    }
}

fn scaleScalarLoop(values: []f64, scalar: f64) void {
    for (values) |*value| {
        value.* *= scalar;
    }
}

fn negateScalarLoop(values: []f64) void {
    for (values) |*value| {
        value.* = -value.*;
    }
}

fn innerProductScalarLoop(lhs: []const f64, rhs: []const f64) f64 {
    var sum: f64 = 0;
    for (lhs, rhs) |left, right| {
        sum += left * right;
    }
    return sum;
}

fn gridVertexIndex(i: u32, j: u32, ny: u32) u32 {
    return i * (ny + 1) + j;
}

fn horizontalEdgeIndex(i: u32, j: u32, nx: u32) u32 {
    return j * nx + i;
}

fn verticalEdgeIndex(i: u32, j: u32, ny: u32, horizontal_edge_count: u32) u32 {
    return horizontal_edge_count + i * ny + j;
}

fn diagonalEdgeIndex(i: u32, j: u32, ny: u32, horizontal_edge_count: u32, vertical_edge_count: u32) u32 {
    return horizontal_edge_count + vertical_edge_count + i * ny + j;
}

fn buildBoundary1Matrix(allocator: std.mem.Allocator, nx: u32, ny: u32) !flux.math.sparse.CsrMatrix(i8) {
    const vertex_count: u32 = (nx + 1) * (ny + 1);
    const horizontal_edge_count: u32 = nx * (ny + 1);
    const vertical_edge_count: u32 = (nx + 1) * ny;
    const diagonal_edge_count: u32 = nx * ny;
    const edge_count: u32 = horizontal_edge_count + vertical_edge_count + diagonal_edge_count;

    var boundary = try flux.math.sparse.CsrMatrix(i8).init(allocator, edge_count, vertex_count, 2 * edge_count);
    errdefer boundary.deinit(allocator);

    var edge_idx: u32 = 0;

    for (0..ny + 1) |j_u| {
        for (0..nx) |i_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            boundary.row_ptr[edge_idx] = 2 * edge_idx;
            boundary.col_idx[2 * edge_idx] = gridVertexIndex(i, j, ny);
            boundary.values[2 * edge_idx] = -1;
            boundary.col_idx[2 * edge_idx + 1] = gridVertexIndex(i + 1, j, ny);
            boundary.values[2 * edge_idx + 1] = 1;
            edge_idx += 1;
        }
    }

    for (0..nx + 1) |i_u| {
        for (0..ny) |j_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            boundary.row_ptr[edge_idx] = 2 * edge_idx;
            boundary.col_idx[2 * edge_idx] = gridVertexIndex(i, j, ny);
            boundary.values[2 * edge_idx] = -1;
            boundary.col_idx[2 * edge_idx + 1] = gridVertexIndex(i, j + 1, ny);
            boundary.values[2 * edge_idx + 1] = 1;
            edge_idx += 1;
        }
    }

    for (0..nx) |i_u| {
        for (0..ny) |j_u| {
            const i: u32 = @intCast(i_u);
            const j: u32 = @intCast(j_u);
            boundary.row_ptr[edge_idx] = 2 * edge_idx;
            boundary.col_idx[2 * edge_idx] = gridVertexIndex(i, j, ny);
            boundary.values[2 * edge_idx] = -1;
            boundary.col_idx[2 * edge_idx + 1] = gridVertexIndex(i + 1, j + 1, ny);
            boundary.values[2 * edge_idx + 1] = 1;
            edge_idx += 1;
        }
    }

    std.debug.assert(edge_idx == edge_count);
    boundary.row_ptr[edge_count] = 2 * edge_count;
    return boundary;
}

fn buildBoundary2Matrix(allocator: std.mem.Allocator, nx: u32, ny: u32) !flux.math.sparse.CsrMatrix(i8) {
    const horizontal_edge_count: u32 = nx * (ny + 1);
    const vertical_edge_count: u32 = (nx + 1) * ny;
    const edge_count: u32 = horizontal_edge_count + vertical_edge_count + nx * ny;
    const face_count: u32 = 2 * nx * ny;

    var boundary = try flux.math.sparse.CsrMatrix(i8).init(allocator, face_count, edge_count, 3 * face_count);
    errdefer boundary.deinit(allocator);

    var face_idx: u32 = 0;
    for (0..nx) |i_u| {
        const i: u32 = @intCast(i_u);
        for (0..ny) |j_u| {
            const j: u32 = @intCast(j_u);
            const h_ij = horizontalEdgeIndex(i, j, nx);
            const h_i_jp1 = horizontalEdgeIndex(i, j + 1, nx);
            const v_ip1_j = verticalEdgeIndex(i + 1, j, ny, horizontal_edge_count);
            const v_i_j = verticalEdgeIndex(i, j, ny, horizontal_edge_count);
            const d_ij = diagonalEdgeIndex(i, j, ny, horizontal_edge_count, vertical_edge_count);

            boundary.row_ptr[face_idx] = 3 * face_idx;
            boundary.col_idx[3 * face_idx + 0] = h_ij;
            boundary.values[3 * face_idx + 0] = 1;
            boundary.col_idx[3 * face_idx + 1] = v_ip1_j;
            boundary.values[3 * face_idx + 1] = 1;
            boundary.col_idx[3 * face_idx + 2] = d_ij;
            boundary.values[3 * face_idx + 2] = -1;
            face_idx += 1;

            boundary.row_ptr[face_idx] = 3 * face_idx;
            boundary.col_idx[3 * face_idx + 0] = h_i_jp1;
            boundary.values[3 * face_idx + 0] = -1;
            boundary.col_idx[3 * face_idx + 1] = v_i_j;
            boundary.values[3 * face_idx + 1] = -1;
            boundary.col_idx[3 * face_idx + 2] = d_ij;
            boundary.values[3 * face_idx + 2] = 1;
            face_idx += 1;
        }
    }

    std.debug.assert(face_idx == face_count);
    boundary.row_ptr[face_count] = 3 * face_count;
    return boundary;
}

fn buildBoundary3LikeMatrix(allocator: std.mem.Allocator, row_count: u32) !flux.math.sparse.CsrMatrix(i8) {
    const col_count: u32 = row_count + 3;
    var boundary = try flux.math.sparse.CsrMatrix(i8).init(allocator, row_count, col_count, 4 * row_count);
    errdefer boundary.deinit(allocator);

    for (0..row_count) |row_idx_usize| {
        const row_idx: u32 = @intCast(row_idx_usize);
        const start = 4 * row_idx;
        boundary.row_ptr[row_idx] = start;
        boundary.col_idx[start + 0] = row_idx;
        boundary.values[start + 0] = 1;
        boundary.col_idx[start + 1] = row_idx + 1;
        boundary.values[start + 1] = -1;
        boundary.col_idx[start + 2] = row_idx + 2;
        boundary.values[start + 2] = 1;
        boundary.col_idx[start + 3] = row_idx + 3;
        boundary.values[start + 3] = -1;
    }
    boundary.row_ptr[row_count] = 4 * row_count;
    return boundary;
}

fn compareIncidenceTransposeMultiply(
    allocator: std.mem.Allocator,
    name: []const u8,
    matrix: flux.math.sparse.CsrMatrix(i8),
    specialized_matrix: flux.math.sparse.PackedIncidenceMatrix,
) !IncidenceComparisonResult {
    var generic_matrix = try flux.math.sparse.PackedIncidenceMatrix.fromCsr(allocator, matrix);
    defer generic_matrix.deinit(allocator);

    const input = try allocator.alloc(f64, matrix.n_rows);
    defer allocator.free(input);
    const dense_output = try allocator.alloc(f64, matrix.n_cols);
    defer allocator.free(dense_output);
    const generic_output = try allocator.alloc(f64, matrix.n_cols);
    defer allocator.free(generic_output);
    const specialized_output = try allocator.alloc(f64, matrix.n_cols);
    defer allocator.free(specialized_output);

    var rng = std.Random.DefaultPrng.init(0x1AC1_D3C3);
    fillRandom(&rng, input);

    for (0..incidence_warmup_iterations) |_| {
        @memset(dense_output, 0.0);
        matrix.transpose_multiply(input, dense_output);
        @memset(generic_output, 0.0);
        generic_matrix.transpose_multiply(input, generic_output);
        @memset(specialized_output, 0.0);
        specialized_matrix.transpose_multiply(input, specialized_output);
    }

    var dense_timings: [incidence_measured_iterations]u64 = undefined;
    var generic_timings: [incidence_measured_iterations]u64 = undefined;
    var specialized_timings: [incidence_measured_iterations]u64 = undefined;
    for (0..incidence_measured_iterations) |sample_idx| {
        const order = switch (sample_idx % 3) {
            0 => [_]u2{ 0, 1, 2 },
            1 => [_]u2{ 1, 2, 0 },
            else => [_]u2{ 2, 0, 1 },
        };
        for (order) |variant_idx| {
            switch (variant_idx) {
                0 => {
                    @memset(dense_output, 0.0);
                    var timer = try std.time.Timer.start();
                    matrix.transpose_multiply(input, dense_output);
                    dense_timings[sample_idx] = timer.read();
                },
                1 => {
                    @memset(generic_output, 0.0);
                    var timer = try std.time.Timer.start();
                    generic_matrix.transpose_multiply(input, generic_output);
                    generic_timings[sample_idx] = timer.read();
                },
                2 => {
                    @memset(specialized_output, 0.0);
                    var timer = try std.time.Timer.start();
                    specialized_matrix.transpose_multiply(input, specialized_output);
                    specialized_timings[sample_idx] = timer.read();
                },
                else => unreachable,
            }
        }
    }

    std.mem.sortUnstable(u64, &dense_timings, {}, std.sort.asc(u64));
    std.mem.sortUnstable(u64, &generic_timings, {}, std.sort.asc(u64));
    std.mem.sortUnstable(u64, &specialized_timings, {}, std.sort.asc(u64));

    const dense_median_ns = dense_timings[incidence_measured_iterations / 2];
    const generic_packed_median_ns = generic_timings[incidence_measured_iterations / 2];
    const specialized_median_ns = specialized_timings[incidence_measured_iterations / 2];

    for (dense_output, generic_output) |expected, actual| {
        try std.testing.expectApproxEqAbs(expected, actual, 1e-12);
    }
    for (dense_output, specialized_output) |expected, actual| {
        try std.testing.expectApproxEqAbs(expected, actual, 1e-12);
    }

    return .{
        .name = name,
        .nnz = matrix.nnz(),
        .dense_value_bytes = matrix.values.len * @sizeOf(i8),
        .generic_packed_value_bytes = generic_matrix.signBytes(),
        .specialized_value_bytes = specialized_matrix.signBytes(),
        .dense_median_ns = dense_median_ns,
        .generic_packed_median_ns = generic_packed_median_ns,
        .specialized_median_ns = specialized_median_ns,
    };
}

fn runIncidenceComparisonBenchmarks(allocator: std.mem.Allocator) ![3]IncidenceComparisonResult {
    var results: [3]IncidenceComparisonResult = undefined;

    var boundary1 = try buildBoundary1Matrix(allocator, incidence_grid_nx, incidence_grid_ny);
    defer boundary1.deinit(allocator);
    var boundary1_specialized = try flux.math.sparse.PackedIncidenceMatrix.fromBoundaryCsr(allocator, 1, boundary1);
    defer boundary1_specialized.deinit(allocator);
    results[0] = try compareIncidenceTransposeMultiply(
        allocator,
        "boundary_1_transpose_multiply",
        boundary1,
        boundary1_specialized,
    );

    var boundary2 = try buildBoundary2Matrix(allocator, incidence_grid_nx, incidence_grid_ny);
    defer boundary2.deinit(allocator);
    var boundary2_specialized = try flux.math.sparse.PackedIncidenceMatrix.fromBoundaryCsr(allocator, 2, boundary2);
    defer boundary2_specialized.deinit(allocator);
    results[1] = try compareIncidenceTransposeMultiply(
        allocator,
        "boundary_2_transpose_multiply",
        boundary2,
        boundary2_specialized,
    );

    const boundary3_rows: u32 = incidence_grid_nx * incidence_grid_ny;
    var boundary3 = try buildBoundary3LikeMatrix(allocator, boundary3_rows);
    defer boundary3.deinit(allocator);
    var boundary3_specialized = try flux.math.sparse.PackedIncidenceMatrix.fromBoundaryCsr(allocator, 3, boundary3);
    defer boundary3_specialized.deinit(allocator);
    results[2] = try compareIncidenceTransposeMultiply(
        allocator,
        "boundary_3_like_transpose_multiply",
        boundary3,
        boundary3_specialized,
    );

    return results;
}

// ── I/O helpers (Zig 0.15 compatible) ───────────────────────────────────

fn stderrWriter() @TypeOf((std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter()) {
    return (std.fs.File{ .handle = std.posix.STDERR_FILENO }).deprecatedWriter();
}

fn stdoutWriter() @TypeOf((std.fs.File{ .handle = std.posix.STDOUT_FILENO }).deprecatedWriter()) {
    return (std.fs.File{ .handle = std.posix.STDOUT_FILENO }).deprecatedWriter();
}

// ── Benchmark implementations ───────────────────────────────────────────

/// d₀: exterior derivative on 0-forms (sparse matrix–vector multiply).
fn benchExteriorDerivativeD0(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.exteriorDerivative(flux.forms.Primal, 0) catch unreachable).apply(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// d₁: exterior derivative on 1-forms.
fn benchExteriorDerivativeD1(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.exteriorDerivative(flux.forms.Primal, 1) catch unreachable).apply(ctx.allocator, ctx.c1) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₀: diagonal Hodge star on 0-forms.
fn benchHodgeStar0(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.hodgeStar(0) catch unreachable).apply(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₁: Whitney mass matrix SpMV on 1-forms.
fn benchHodgeStar1(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.hodgeStar(1) catch unreachable).apply(ctx.allocator, ctx.c1) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₂: diagonal Hodge star on 2-forms.
fn benchHodgeStar2(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.hodgeStar(2) catch unreachable).apply(ctx.allocator, ctx.c2) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★⁻¹₁: inverse Whitney Hodge star via preconditioned CG solve.
fn benchHodgeStarInverse1(ctx: *BenchmarkContext) void {
    const DualC1 = flux.forms.Cochain(Mesh2D, 1, flux.forms.Dual);
    var dual = DualC1.init(ctx.allocator, ctx.mesh) catch unreachable;
    defer dual.deinit(ctx.allocator);
    var rng = std.Random.DefaultPrng.init(0xDEAD);
    fillRandom(&rng, dual.values);

    var result = (ctx.operator_context.hodgeStarInverse(1) catch unreachable).apply(ctx.allocator, dual) catch unreachable;
    result.deinit(ctx.allocator);
}

/// Δ₀: repeated application via a reused assembled operator.
fn benchLaplacian0(ctx: *BenchmarkContext) void {
    var result = (ctx.operator_context.laplacian(0) catch unreachable).apply(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// Δ₀: reference composed path without stored assembly.
fn benchLaplacian0Composed(ctx: *BenchmarkContext) void {
    var result = flux.operators.laplacian.laplacian_composed(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// Cochain add: pointwise addition of two 1-cochains.
fn benchCochainAdd(ctx: *BenchmarkContext) void {
    ctx.c1.add(ctx.c1_other);
}

/// Cochain scale: scalar multiplication on a 1-cochain.
fn benchCochainScale(ctx: *BenchmarkContext) void {
    ctx.c1.scale(arithmetic_scale_factor);
}

/// Cochain negate: sign flip on a 1-cochain.
fn benchCochainNegate(ctx: *BenchmarkContext) void {
    ctx.c1.negate();
}

/// Cochain inner product: ⟨a, b⟩ on 0-cochains.
fn benchCochainInnerProduct(ctx: *BenchmarkContext) void {
    const result = ctx.c0.inner_product(ctx.c0_other);
    std.mem.doNotOptimizeAway(result);
}

fn benchMaxwellCavityStep256(ctx: *BenchmarkContext) void {
    maxwell.step(2, ctx.allocator, &ctx.cavity_state, cavityDt()) catch unreachable;
}

fn benchSurfaceLaplacianAssemblyDirect(ctx: *BenchmarkContext) void {
    const operator_context = SurfaceOperatorContext.init(ctx.allocator, ctx.surface_mesh) catch unreachable;
    defer operator_context.deinit();
    _ = operator_context.laplacian(0) catch unreachable;
}

fn benchDiffusionSurfaceSystemInit(ctx: *BenchmarkContext) void {
    const dt = diffusion_surface_final_time / @as(f64, @floatFromInt(diffusion_surface_steps));
    var mesh = SurfaceMesh.sphere(ctx.allocator, 1.0, diffusion_surface_refinement) catch unreachable;
    defer mesh.deinit(ctx.allocator);

    var system = diffusion_sphere.SurfaceSystem.init(ctx.allocator, &mesh, dt) catch unreachable;
    defer system.deinit(ctx.allocator);
}

fn arithmeticCompareDefs() [arithmetic_case_lengths.len * 4]ArithmeticCompareDef {
    var defs: [arithmetic_case_lengths.len * 4]ArithmeticCompareDef = undefined;
    var i = 0;

    inline for (arithmetic_case_labels, 0..) |label, case_index| {
        defs[i] = .{
            .name = std.fmt.comptimePrint("cochain_add_{s}", .{label}),
            .operation = .add,
            .case_index = case_index,
        };
        i += 1;
        defs[i] = .{
            .name = std.fmt.comptimePrint("cochain_scale_{s}", .{label}),
            .operation = .scale,
            .case_index = case_index,
        };
        i += 1;
        defs[i] = .{
            .name = std.fmt.comptimePrint("cochain_negate_{s}", .{label}),
            .operation = .negate,
            .case_index = case_index,
        };
        i += 1;
        defs[i] = .{
            .name = std.fmt.comptimePrint("cochain_inner_product_{s}", .{label}),
            .operation = .inner_product,
            .case_index = case_index,
        };
        i += 1;
    }

    return defs;
}

const arithmetic_compare_defs = arithmeticCompareDefs();

const default_benchmarks = [_]BenchmarkDef{
    .{ .name = "exterior_derivative_d0", .run = benchExteriorDerivativeD0, .scope = .public_library },
    .{ .name = "exterior_derivative_d1", .run = benchExteriorDerivativeD1, .scope = .public_library },
    .{ .name = "hodge_star_0", .run = benchHodgeStar0, .scope = .public_library },
    .{ .name = "hodge_star_1_whitney", .run = benchHodgeStar1, .scope = .public_library },
    .{ .name = "hodge_star_2", .run = benchHodgeStar2, .scope = .public_library },
    .{ .name = "hodge_star_inverse_1_cg", .run = benchHodgeStarInverse1, .scope = .public_library },
    .{ .name = "laplacian_0", .run = benchLaplacian0, .scope = .public_library },
    .{ .name = "laplacian_0_composed", .run = benchLaplacian0Composed, .scope = .public_library },
    .{
        .name = "cochain_add",
        .run = benchCochainAdd,
        .minimum_repetitions = small_cochain_repetitions,
        .version = micro_benchmark_method_version,
        .class = .info,
        .scope = .public_library,
    },
    .{
        .name = "cochain_scale",
        .run = benchCochainScale,
        .minimum_repetitions = small_cochain_repetitions,
        .version = micro_benchmark_method_version,
        .class = .info,
        .scope = .public_library,
    },
    .{
        .name = "cochain_negate",
        .run = benchCochainNegate,
        .minimum_repetitions = small_cochain_repetitions,
        .version = micro_benchmark_method_version,
        .class = .info,
        .scope = .public_library,
    },
    .{
        .name = "cochain_inner_product",
        .run = benchCochainInnerProduct,
        .minimum_repetitions = small_cochain_repetitions,
        .version = micro_benchmark_method_version,
        .class = .info,
        .scope = .public_library,
    },
    .{ .name = "maxwell_cavity_step_256", .run = benchMaxwellCavityStep256, .scope = .shipped_example },
    .{
        .name = "surface_laplacian_0_assembly_direct",
        .run = benchSurfaceLaplacianAssemblyDirect,
        .class = .info,
        .scope = .public_library,
    },
    .{
        .name = "diffusion_surface_system_init",
        .run = benchDiffusionSurfaceSystemInit,
        .class = .info,
        .scope = .shipped_example,
    },
};

// ── Runner ──────────────────────────────────────────────────────────────

fn runBenchmark(def: BenchmarkDef, ctx: *BenchmarkContext) BenchmarkResult {
    // Warmup — let branch predictors and caches settle.
    for (0..warmup_iterations) |_| {
        for (0..def.minimum_repetitions) |_| {
            def.run(ctx);
        }
    }

    const repetitions = calibrateRepetitions(def, ctx);

    // Measured iterations.
    var timings: [measured_iterations]u64 = undefined;
    for (&timings) |*t| {
        var timer = std.time.Timer.start() catch unreachable;
        for (0..repetitions) |_| {
            def.run(ctx);
        }
        t.* = @divFloor(timer.read(), repetitions);
    }

    // Sort for median/min/max.
    std.mem.sortUnstable(u64, &timings, {}, std.sort.asc(u64));

    return .{
        .name = def.name,
        .version = def.version,
        .class = def.class,
        .repetitions = repetitions,
        .median_ns = timings[measured_iterations / 2],
        .min_ns = timings[0],
        .max_ns = timings[measured_iterations - 1],
        .iterations = measured_iterations,
    };
}

fn calibrateRepetitions(def: BenchmarkDef, ctx: *BenchmarkContext) u32 {
    var repetitions = def.minimum_repetitions;
    while (true) {
        var timer = std.time.Timer.start() catch unreachable;
        for (0..repetitions) |_| {
            def.run(ctx);
        }
        const elapsed_ns = timer.read();
        if (elapsed_ns >= target_sample_duration_ns) return repetitions;
        if (repetitions >= max_sample_repetitions) return repetitions;

        const doubled = @as(u64, repetitions) * 2;
        repetitions = @intCast(@min(doubled, max_sample_repetitions));
    }
}

fn operationLabel(operation: ArithmeticOperation) []const u8 {
    return switch (operation) {
        .add => "add",
        .scale => "scale",
        .negate => "negate",
        .inner_product => "inner_product",
    };
}

const ArithmeticVariant = enum {
    scalar,
    default_impl,
};

fn resetArithmeticCase(arithmetic_case: *ArithmeticCase, lhs_seed: []const f64, rhs_seed: []const f64) void {
    @memcpy(arithmetic_case.lhs.values, lhs_seed);
    @memcpy(arithmetic_case.rhs.values, rhs_seed);
}

fn runArithmeticVariant(
    arithmetic_case: *ArithmeticCase,
    operation: ArithmeticOperation,
    variant: ArithmeticVariant,
) void {
    switch (operation) {
        .add => switch (variant) {
            .scalar => addScalarLoop(arithmetic_case.lhs.values, arithmetic_case.rhs.values),
            .default_impl => arithmetic_case.lhs.add(arithmetic_case.rhs),
        },
        .scale => switch (variant) {
            .scalar => scaleScalarLoop(arithmetic_case.lhs.values, arithmetic_scale_factor),
            .default_impl => arithmetic_case.lhs.scale(arithmetic_scale_factor),
        },
        .negate => switch (variant) {
            .scalar => negateScalarLoop(arithmetic_case.lhs.values),
            .default_impl => arithmetic_case.lhs.negate(),
        },
        .inner_product => {
            const result = switch (variant) {
                .scalar => innerProductScalarLoop(arithmetic_case.lhs.values, arithmetic_case.rhs.values),
                .default_impl => arithmetic_case.lhs.inner_product(arithmetic_case.rhs),
            };
            std.mem.doNotOptimizeAway(result);
        },
    }
}

fn measureArithmeticComparison(def: ArithmeticCompareDef, ctx: *BenchmarkContext) !ArithmeticComparisonResult {
    const arithmetic_case = &ctx.arithmetic_cases[def.case_index];
    const lhs_seed = try ctx.allocator.dupe(f64, arithmetic_case.lhs.values);
    defer ctx.allocator.free(lhs_seed);
    const rhs_seed = try ctx.allocator.dupe(f64, arithmetic_case.rhs.values);
    defer ctx.allocator.free(rhs_seed);

    for (0..warmup_iterations) |warmup_idx| {
        const variants = if (warmup_idx % 2 == 0)
            [_]ArithmeticVariant{ .scalar, .default_impl }
        else
            [_]ArithmeticVariant{ .default_impl, .scalar };
        for (variants) |variant| {
            resetArithmeticCase(arithmetic_case, lhs_seed, rhs_seed);
            runArithmeticVariant(arithmetic_case, def.operation, variant);
        }
    }

    var scalar_timings: [measured_iterations]u64 = undefined;
    var default_timings: [measured_iterations]u64 = undefined;
    for (0..measured_iterations) |sample_idx| {
        const variants = if (sample_idx % 2 == 0)
            [_]ArithmeticVariant{ .scalar, .default_impl }
        else
            [_]ArithmeticVariant{ .default_impl, .scalar };
        for (variants) |variant| {
            resetArithmeticCase(arithmetic_case, lhs_seed, rhs_seed);
            var timer = try std.time.Timer.start();
            runArithmeticVariant(arithmetic_case, def.operation, variant);
            const elapsed_ns = timer.read();
            switch (variant) {
                .scalar => scalar_timings[sample_idx] = elapsed_ns,
                .default_impl => default_timings[sample_idx] = elapsed_ns,
            }
        }
    }

    std.mem.sortUnstable(u64, &scalar_timings, {}, std.sort.asc(u64));
    std.mem.sortUnstable(u64, &default_timings, {}, std.sort.asc(u64));

    return .{
        .name = def.name,
        .operation = def.operation,
        .case_label = arithmetic_case_labels[def.case_index],
        .version = def.version,
        .scalar_median_ns = scalar_timings[measured_iterations / 2],
        .default_median_ns = default_timings[measured_iterations / 2],
    };
}

fn runArithmeticCompareBenchmarks(ctx: *BenchmarkContext) ![arithmetic_compare_defs.len]ArithmeticComparisonResult {
    var results: [arithmetic_compare_defs.len]ArithmeticComparisonResult = undefined;
    for (arithmetic_compare_defs, 0..) |def, i| {
        results[i] = try measureArithmeticComparison(def, ctx);
    }
    return results;
}

// ── Baseline I/O ────────────────────────────────────────────────────────

const Baseline = struct {
    name: []const u8,
    version: u32 = 1,
    class: BenchmarkClass = .gate,
    median_ns: u64,
};

const BaselineFile = struct {
    suite_version: u32 = 1,
    mesh_size: []const u8,
    iterations: u32,
    benchmarks: []const Baseline,
};

fn readBaselines(allocator: std.mem.Allocator) !?std.json.Parsed(BaselineFile) {
    const file = std.fs.cwd().openFile("bench/baselines.json", .{}) catch |err| {
        if (err == error.FileNotFound) return null;
        return err;
    };
    defer file.close();

    const content = try file.readToEndAlloc(allocator, 1024 * 1024);
    defer allocator.free(content);

    return try std.json.parseFromSlice(BaselineFile, allocator, content, .{
        .allocate = .alloc_always,
    });
}

fn writeBaselines(results: []const BenchmarkResult) !void {
    var file = try std.fs.cwd().createFile("bench/baselines.json", .{});
    defer file.close();

    const writer = file.deprecatedWriter();
    try writer.writeAll("{\n");
    try writer.print("  \"suite_version\": {d},\n", .{benchmark_suite_version});
    try writer.print("  \"mesh_size\": \"{d}x{d}\",\n", .{ grid_nx, grid_ny });
    try writer.print("  \"iterations\": {d},\n", .{measured_iterations});
    try writer.writeAll("  \"benchmarks\": [\n");

    for (results, 0..) |r, i| {
        try writer.print(
            "    {{\"name\": \"{s}\", \"version\": {d}, \"class\": \"{s}\", \"median_ns\": {d}}}",
            .{ r.name, r.version, @tagName(r.class), r.median_ns },
        );
        if (i < results.len - 1) {
            try writer.writeAll(",\n");
        } else {
            try writer.writeAll("\n");
        }
    }

    try writer.writeAll("  ]\n");
    try writer.writeAll("}\n");
}

// ── Output formatting ───────────────────────────────────────────────────

fn printTable(
    writer: anytype,
    results: []const BenchmarkResult,
    baselines: ?[]const Baseline,
) !void {
    try writer.writeAll("\n");
    try writer.print("  flux benchmark suite — {d}x{d} mesh, {d} iterations\n", .{
        grid_nx, grid_ny, measured_iterations,
    });
    try writer.print("  suite version: {d}\n", .{benchmark_suite_version});
    try writer.writeAll("  ─────────────────────────────────────────────────────────────────\n");

    if (baselines != null) {
        try writer.print("  {s:<30} {s:>12} {s:>12} {s:>12} {s:>8}\n", .{
            "benchmark", "median", "min", "max", "vs base",
        });
    } else {
        try writer.print("  {s:<30} {s:>12} {s:>12} {s:>12}\n", .{
            "benchmark", "median", "min", "max",
        });
    }

    try writer.writeAll("  ─────────────────────────────────────────────────────────────────\n");

    for (results) |r| {
        const baseline = if (baselines) |bl| findBaseline(bl, r.name) else null;

        try writer.print("  {s:<30} ", .{r.name});
        try printDuration(writer, r.median_ns);
        try writer.writeAll("  ");
        try printDuration(writer, r.min_ns);
        try writer.writeAll("  ");
        try printDuration(writer, r.max_ns);

        if (baseline) |base| {
            if (base.version != r.version) {
                try writer.print("  method v{d}→v{d}", .{ base.version, r.version });
            } else if (r.class == .info) {
                try writer.writeAll("  info only");
            } else {
                const ratio = @as(f64, @floatFromInt(r.median_ns)) / @as(f64, @floatFromInt(base.median_ns));
                const pct = (ratio - 1.0) * 100.0;
                if (pct > regression_threshold * 100.0) {
                    try writer.print("  +{d:.1}% REGRESSED", .{pct});
                } else if (pct < -5.0) {
                    try writer.print("  {d:.1}% faster", .{-pct});
                } else if (pct >= 0) {
                    try writer.print("  +{d:.1}%", .{pct});
                } else {
                    try writer.print("  {d:.1}%", .{pct});
                }
            }
        }

        try writer.writeAll("\n");
    }

    try writer.writeAll("  ─────────────────────────────────────────────────────────────────\n\n");
}

fn printArithmeticComparisons(writer: anytype, results: []const ArithmeticComparisonResult) !void {
    try writer.writeAll("  Arithmetic same-run comparisons (default cochain path vs scalar loop)\n");
    try writer.writeAll("  ───────────────────────────────────────────────────────────────────────────────────────\n");
    try writer.print("  {s:<16} {s:<8} {s:>12} {s:>12} {s:>10}\n", .{
        "operation",
        "case",
        "scalar",
        "default",
        "speedup",
    });
    try writer.writeAll("  ───────────────────────────────────────────────────────────────────────────────────────\n");

    for (results) |result| {
        const speedup = @as(f64, @floatFromInt(result.scalar_median_ns)) /
            @as(f64, @floatFromInt(result.default_median_ns));
        try writer.print("  {s:<16} {s:<8} ", .{
            operationLabel(result.operation),
            result.case_label,
        });
        try printDuration(writer, result.scalar_median_ns);
        try writer.writeAll("  ");
        try printDuration(writer, result.default_median_ns);
        try writer.print("  {d:>7.2}x\n", .{speedup});
    }

    try writer.writeAll("  ───────────────────────────────────────────────────────────────────────────────────────\n\n");
}

fn printIncidenceComparisons(writer: anytype, results: []const IncidenceComparisonResult) !void {
    try writer.print("  Incidence transpose-multiply comparison ({d}x{d} grid)\n", .{
        incidence_grid_nx,
        incidence_grid_ny,
    });
    try writer.writeAll("  ──────────────────────────────────────────────────────────────────────────────────────────────────\n");
    try writer.print("  {s:<32} {s:>10} {s:>10} {s:>10} {s:>8} {s:>8} {s:>8} {s:>8}\n", .{
        "benchmark",
        "i8",
        "generic",
        "special",
        "gen x",
        "spec x",
        "gen b/e",
        "spec b/e",
    });
    try writer.writeAll("  ──────────────────────────────────────────────────────────────────────────────────────────────────\n");

    for (results) |result| {
        const generic_speedup = @as(f64, @floatFromInt(result.dense_median_ns)) /
            @as(f64, @floatFromInt(result.generic_packed_median_ns));
        const specialized_speedup = @as(f64, @floatFromInt(result.dense_median_ns)) /
            @as(f64, @floatFromInt(result.specialized_median_ns));
        const generic_bits_per_entry = 8.0 *
            @as(f64, @floatFromInt(result.generic_packed_value_bytes)) /
            @as(f64, @floatFromInt(result.nnz));
        const specialized_bits_per_entry = 8.0 *
            @as(f64, @floatFromInt(result.specialized_value_bytes)) /
            @as(f64, @floatFromInt(result.nnz));

        try writer.print("  {s:<32} ", .{result.name});
        try printDuration(writer, result.dense_median_ns);
        try writer.writeAll("  ");
        try printDuration(writer, result.generic_packed_median_ns);
        try writer.writeAll("  ");
        try printDuration(writer, result.specialized_median_ns);
        try writer.print("  {d:>6.2}x  {d:>7.2}x  {d:>8.3}  {d:>8.3}\n", .{
            generic_speedup,
            specialized_speedup,
            generic_bits_per_entry,
            specialized_bits_per_entry,
        });
    }

    try writer.writeAll("  ──────────────────────────────────────────────────────────────────────────────────────────────────\n\n");
}

fn printDuration(writer: anytype, ns: u64) !void {
    if (ns < 1_000) {
        try writer.print("{d:>8} ns", .{ns});
    } else if (ns < 1_000_000) {
        const us = @as(f64, @floatFromInt(ns)) / 1_000.0;
        try writer.print("{d:>8.1} µs", .{us});
    } else if (ns < 1_000_000_000) {
        const ms = @as(f64, @floatFromInt(ns)) / 1_000_000.0;
        try writer.print("{d:>8.2} ms", .{ms});
    } else {
        const s = @as(f64, @floatFromInt(ns)) / 1_000_000_000.0;
        try writer.print("{d:>8.3}  s", .{s});
    }
}

fn printJson(writer: anytype, results: []const BenchmarkResult) !void {
    try writer.writeAll("{\n");
    try writer.print("  \"suite_version\": {d},\n", .{benchmark_suite_version});
    try writer.writeAll("  \"mode\": \"default\",\n");
    try writer.print("  \"mesh_size\": \"{d}x{d}\",\n", .{ grid_nx, grid_ny });
    try writer.print("  \"iterations\": {d},\n", .{measured_iterations});
    try writer.writeAll("  \"context\": {\n");
    try writer.print("    \"zig_version\": \"{s}\",\n", .{builtin.zig_version_string});
    try writer.print("    \"arch\": \"{s}\",\n", .{@tagName(builtin.target.cpu.arch)});
    try writer.print("    \"os\": \"{s}\",\n", .{@tagName(builtin.target.os.tag)});
    try writer.print("    \"target_sample_duration_ns\": {d}\n", .{target_sample_duration_ns});
    try writer.writeAll("  },\n");
    try writer.writeAll("  \"benchmarks\": [\n");

    for (results, 0..) |r, i| {
        try writer.print(
            "    {{\"name\": \"{s}\", \"version\": {d}, \"class\": \"{s}\", \"repetitions\": {d}, \"median_ns\": {d}, \"min_ns\": {d}, \"max_ns\": {d}}}",
            .{ r.name, r.version, @tagName(r.class), r.repetitions, r.median_ns, r.min_ns, r.max_ns },
        );
        if (i < results.len - 1) {
            try writer.writeAll(",\n");
        } else {
            try writer.writeAll("\n");
        }
    }

    try writer.writeAll("  ]\n");
    try writer.writeAll("}\n");
}

fn printCompareJson(
    writer: anytype,
    arithmetic_results: []const ArithmeticComparisonResult,
    incidence_results: []const IncidenceComparisonResult,
) !void {
    try writer.writeAll("{\n");
    try writer.print("  \"suite_version\": {d},\n", .{benchmark_suite_version});
    try writer.writeAll("  \"mode\": \"compare\",\n");
    try writer.writeAll("  \"context\": {\n");
    try writer.print("    \"zig_version\": \"{s}\",\n", .{builtin.zig_version_string});
    try writer.print("    \"arch\": \"{s}\",\n", .{@tagName(builtin.target.cpu.arch)});
    try writer.print("    \"os\": \"{s}\",\n", .{@tagName(builtin.target.os.tag)});
    try writer.print("    \"measured_iterations\": {d}\n", .{measured_iterations});
    try writer.writeAll("  },\n");
    try writer.writeAll("  \"arithmetic_comparisons\": [\n");
    for (arithmetic_results, 0..) |result, i| {
        try writer.print(
            "    {{\"name\": \"{s}\", \"operation\": \"{s}\", \"case\": \"{s}\", \"version\": {d}, \"scalar_median_ns\": {d}, \"default_median_ns\": {d}}}",
            .{
                result.name,
                operationLabel(result.operation),
                result.case_label,
                result.version,
                result.scalar_median_ns,
                result.default_median_ns,
            },
        );
        if (i < arithmetic_results.len - 1) {
            try writer.writeAll(",\n");
        } else {
            try writer.writeAll("\n");
        }
    }
    try writer.writeAll("  ],\n");
    try writer.writeAll("  \"incidence_comparisons\": [\n");
    for (incidence_results, 0..) |result, i| {
        try writer.print(
            "    {{\"name\": \"{s}\", \"version\": {d}, \"nnz\": {d}, \"dense_median_ns\": {d}, \"generic_packed_median_ns\": {d}, \"specialized_median_ns\": {d}, \"dense_value_bytes\": {d}, \"generic_packed_value_bytes\": {d}, \"specialized_value_bytes\": {d}}}",
            .{
                result.name,
                result.version,
                result.nnz,
                result.dense_median_ns,
                result.generic_packed_median_ns,
                result.specialized_median_ns,
                result.dense_value_bytes,
                result.generic_packed_value_bytes,
                result.specialized_value_bytes,
            },
        );
        if (i < incidence_results.len - 1) {
            try writer.writeAll(",\n");
        } else {
            try writer.writeAll("\n");
        }
    }
    try writer.writeAll("  ]\n");
    try writer.writeAll("}\n");
}

fn findBaseline(baselines: []const Baseline, name: []const u8) ?Baseline {
    for (baselines) |b| {
        if (std.mem.eql(u8, b.name, name)) return b;
    }
    return null;
}

fn checkRegressions(results: []const BenchmarkResult, baselines: []const Baseline) bool {
    var any_regression = false;
    const stderr = stderrWriter();

    for (results) |r| {
        if (r.class != .gate) continue;
        const base = findBaseline(baselines, r.name) orelse continue;
        if (base.version != r.version) continue;
        if (base.class != .gate) continue;
        const ratio = @as(f64, @floatFromInt(r.median_ns)) / @as(f64, @floatFromInt(base.median_ns));
        if (ratio > 1.0 + regression_threshold) {
            stderr.print("  REGRESSION: {s} — {d:.1}% slower (baseline: ", .{
                r.name, (ratio - 1.0) * 100.0,
            }) catch {};
            printDuration(stderr, base.median_ns) catch {};
            stderr.writeAll(", current: ") catch {};
            printDuration(stderr, r.median_ns) catch {};
            stderr.writeAll(")\n") catch {};
            any_regression = true;
        }
    }

    return any_regression;
}

fn parseBenchmarkArgs(args: []const []const u8) !BenchmarkArgs {
    var parsed: BenchmarkArgs = .{};
    for (args) |arg| {
        if (std.mem.eql(u8, arg, "--check")) {
            parsed.check_mode = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--update")) {
            parsed.update_mode = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--json")) {
            parsed.json_mode = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--compare")) {
            parsed.mode = .compare;
            continue;
        }
        return error.InvalidArguments;
    }

    if (parsed.mode == .compare and (parsed.check_mode or parsed.update_mode)) {
        return error.InvalidArguments;
    }

    return parsed;
}

test "default benchmark suite excludes compare experiments" {
    for (default_benchmarks) |def| {
        try std.testing.expect(def.scope != .compare_experiment);
    }
}

test "compare definitions stay out of the default suite" {
    for (arithmetic_compare_defs) |def| {
        try std.testing.expect(def.scope == .compare_experiment);
    }
}

test "compare mode rejects baseline operations" {
    try std.testing.expectError(error.InvalidArguments, parseBenchmarkArgs(&.{ "--compare", "--check" }));
    try std.testing.expectError(error.InvalidArguments, parseBenchmarkArgs(&.{ "--compare", "--update" }));
}

// ── Main ────────────────────────────────────────────────────────────────

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);
    const benchmark_args = parseBenchmarkArgs(args[1..]) catch {
        try stderrWriter().writeAll(
            "  usage: zig build bench -- [--check] [--update] [--json] [--compare]\n" ++
                "  note: --compare cannot be combined with --check or --update\n",
        );
        return error.InvalidArguments;
    };

    // Build mesh.
    const stderr = stderrWriter();
    try stderr.print("  Building {d}x{d} mesh...\n", .{ grid_nx, grid_ny });

    var mesh = try Mesh2D.uniform_grid(allocator, grid_nx, grid_ny, grid_width, grid_height);
    defer mesh.deinit(allocator);

    try stderr.print("  Mesh: {d} vertices, {d} edges, {d} faces\n", .{
        mesh.num_vertices(), mesh.num_edges(), mesh.num_faces(),
    });

    // Initialize benchmark context with pre-allocated cochains.
    var ctx = try BenchmarkContext.init(allocator, &mesh);
    defer ctx.deinit();

    switch (benchmark_args.mode) {
        .default => {
            var results: [default_benchmarks.len]BenchmarkResult = undefined;
            for (default_benchmarks, 0..) |def, i| {
                try stderr.print("  Running {s}...\n", .{def.name});
                results[i] = runBenchmark(def, &ctx);
            }

            var parsed_baselines = try readBaselines(allocator);
            defer if (parsed_baselines) |*p| p.deinit();

            const baseline_list: ?[]const Baseline = if (parsed_baselines) |p| p.value.benchmarks else null;

            if (benchmark_args.json_mode) {
                try printJson(stdoutWriter(), &results);
            }
            try printTable(stderr, &results, baseline_list);

            if (benchmark_args.update_mode) {
                try writeBaselines(&results);
                try stderr.writeAll("  Baselines updated in bench/baselines.json\n\n");
            }

            if (benchmark_args.check_mode) {
                if (baseline_list) |bl| {
                    if (checkRegressions(&results, bl)) {
                        try stderr.print("  FAIL: one or more benchmarks regressed >{d:.0}%\n\n", .{regression_threshold * 100.0});
                        std.process.exit(1);
                    } else {
                        try stderr.writeAll("  PASS: no regressions detected\n\n");
                    }
                } else {
                    try stderr.writeAll("  SKIP: no baselines.json found, nothing to compare against\n\n");
                }
            }
        },
        .compare => {
            for (arithmetic_compare_defs) |def| {
                try stderr.print("  Running {s} compare...\n", .{def.name});
            }
            const arithmetic_results = try runArithmeticCompareBenchmarks(&ctx);

            try stderr.print("  Building incidence comparison cases on a {d}x{d} grid...\n", .{
                incidence_grid_nx,
                incidence_grid_ny,
            });
            const incidence_results = try runIncidenceComparisonBenchmarks(allocator);

            if (benchmark_args.json_mode) {
                try printCompareJson(stdoutWriter(), &arithmetic_results, &incidence_results);
            }
            try printArithmeticComparisons(stderr, &arithmetic_results);
            try printIncidenceComparisons(stderr, &incidence_results);
        },
    }
}
