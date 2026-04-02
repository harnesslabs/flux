//! Benchmark suite for flux operators.
//!
//! Measures wall-clock performance of key operators on a moderately-sized mesh.
//! Reports median timings (robust to outliers) in both human-readable table
//! and machine-readable JSON formats.
//!
//! Usage:
//!   zig build bench                    — run benchmarks, print results
//!   zig build bench -- --check         — compare against baselines, exit 1 on regression
//!   zig build bench -- --update        — run benchmarks and overwrite baselines.json
//!   zig build bench -- --json          — print JSON to stdout (default: table to stderr)

const std = @import("std");
const flux = @import("flux");

const Mesh2D = flux.Mesh(2, 2);
const PrimalC0 = flux.Cochain(Mesh2D, 0, flux.Primal);
const PrimalC1 = flux.Cochain(Mesh2D, 1, flux.Primal);
const PrimalC2 = flux.Cochain(Mesh2D, 2, flux.Primal);

// ── Configuration ───────────────────────────────────────────────────────

const grid_nx: u32 = 50;
const grid_ny: u32 = 50;
const grid_width: f64 = 50.0;
const grid_height: f64 = 50.0;
const warmup_iterations: u32 = 10;
const measured_iterations: u32 = 100;
/// Maximum allowed regression before CI fails (0.20 = 20%).
const regression_threshold: f64 = 0.20;

// ── Benchmark definition ────────────────────────────────────────────────

const BenchmarkFn = *const fn (*BenchmarkContext) void;

const BenchmarkDef = struct {
    name: []const u8,
    run: BenchmarkFn,
};

const BenchmarkResult = struct {
    name: []const u8,
    median_ns: u64,
    min_ns: u64,
    max_ns: u64,
    iterations: u32,
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
        };
    }

    fn deinit(self: *BenchmarkContext) void {
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
    var result = flux.exterior_derivative(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// d₁: exterior derivative on 1-forms.
fn benchExteriorDerivativeD1(ctx: *BenchmarkContext) void {
    var result = flux.exterior_derivative(ctx.allocator, ctx.c1) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₀: diagonal Hodge star on 0-forms.
fn benchHodgeStar0(ctx: *BenchmarkContext) void {
    var result = flux.hodge_star(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₁: Whitney mass matrix SpMV on 1-forms.
fn benchHodgeStar1(ctx: *BenchmarkContext) void {
    var result = flux.hodge_star(ctx.allocator, ctx.c1) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★₂: diagonal Hodge star on 2-forms.
fn benchHodgeStar2(ctx: *BenchmarkContext) void {
    var result = flux.hodge_star(ctx.allocator, ctx.c2) catch unreachable;
    result.deinit(ctx.allocator);
}

/// ★⁻¹₁: inverse Whitney Hodge star via preconditioned CG solve.
fn benchHodgeStarInverse1(ctx: *BenchmarkContext) void {
    const DualC1 = flux.Cochain(Mesh2D, 1, flux.Dual);
    var dual = DualC1.init(ctx.allocator, ctx.mesh) catch unreachable;
    defer dual.deinit(ctx.allocator);
    var rng = std.Random.DefaultPrng.init(0xDEAD);
    fillRandom(&rng, dual.values);

    var result = flux.hodge_star_inverse(ctx.allocator, dual) catch unreachable;
    result.deinit(ctx.allocator);
}

/// Δ₀: full Laplacian on 0-forms (d, ★, transpose multiply, ★⁻¹ pipeline).
fn benchLaplacian0(ctx: *BenchmarkContext) void {
    var result = flux.laplacian(ctx.allocator, ctx.c0) catch unreachable;
    result.deinit(ctx.allocator);
}

/// Cochain add: pointwise addition of two 1-cochains.
fn benchCochainAdd(ctx: *BenchmarkContext) void {
    for (ctx.c1.values, ctx.c1_other.values) |*a, b| {
        a.* += b;
    }
}

/// Cochain scale: scalar multiplication on a 1-cochain.
fn benchCochainScale(ctx: *BenchmarkContext) void {
    for (ctx.c1.values) |*v| {
        v.* *= 2.5;
    }
}

/// Cochain inner product: ⟨a, b⟩ on 0-cochains.
fn benchCochainInnerProduct(ctx: *BenchmarkContext) void {
    const result = ctx.c0.inner_product(ctx.c0_other);
    std.mem.doNotOptimizeAway(result);
}

const all_benchmarks = [_]BenchmarkDef{
    .{ .name = "exterior_derivative_d0", .run = benchExteriorDerivativeD0 },
    .{ .name = "exterior_derivative_d1", .run = benchExteriorDerivativeD1 },
    .{ .name = "hodge_star_0", .run = benchHodgeStar0 },
    .{ .name = "hodge_star_1_whitney", .run = benchHodgeStar1 },
    .{ .name = "hodge_star_2", .run = benchHodgeStar2 },
    .{ .name = "hodge_star_inverse_1_cg", .run = benchHodgeStarInverse1 },
    .{ .name = "laplacian_0", .run = benchLaplacian0 },
    .{ .name = "cochain_add", .run = benchCochainAdd },
    .{ .name = "cochain_scale", .run = benchCochainScale },
    .{ .name = "cochain_inner_product", .run = benchCochainInnerProduct },
};

// ── Runner ──────────────────────────────────────────────────────────────

fn runBenchmark(def: BenchmarkDef, ctx: *BenchmarkContext) BenchmarkResult {
    // Warmup — let branch predictors and caches settle.
    for (0..warmup_iterations) |_| {
        def.run(ctx);
    }

    // Measured iterations.
    var timings: [measured_iterations]u64 = undefined;
    for (&timings) |*t| {
        var timer = std.time.Timer.start() catch unreachable;
        def.run(ctx);
        t.* = timer.read();
    }

    // Sort for median/min/max.
    std.mem.sortUnstable(u64, &timings, {}, std.sort.asc(u64));

    return .{
        .name = def.name,
        .median_ns = timings[measured_iterations / 2],
        .min_ns = timings[0],
        .max_ns = timings[measured_iterations - 1],
        .iterations = measured_iterations,
    };
}

// ── Baseline I/O ────────────────────────────────────────────────────────

const Baseline = struct {
    name: []const u8,
    median_ns: u64,
};

const BaselineFile = struct {
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
    try writer.print("  \"mesh_size\": \"{d}x{d}\",\n", .{ grid_nx, grid_ny });
    try writer.print("  \"iterations\": {d},\n", .{measured_iterations});
    try writer.writeAll("  \"benchmarks\": [\n");

    for (results, 0..) |r, i| {
        try writer.print("    {{\"name\": \"{s}\", \"median_ns\": {d}}}", .{ r.name, r.median_ns });
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
        const baseline_ns = if (baselines) |bl| findBaseline(bl, r.name) else null;

        try writer.print("  {s:<30} ", .{r.name});
        try printDuration(writer, r.median_ns);
        try writer.writeAll("  ");
        try printDuration(writer, r.min_ns);
        try writer.writeAll("  ");
        try printDuration(writer, r.max_ns);

        if (baseline_ns) |base| {
            const ratio = @as(f64, @floatFromInt(r.median_ns)) / @as(f64, @floatFromInt(base));
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

        try writer.writeAll("\n");
    }

    try writer.writeAll("  ─────────────────────────────────────────────────────────────────\n\n");
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
    try writer.print("  \"mesh_size\": \"{d}x{d}\",\n", .{ grid_nx, grid_ny });
    try writer.print("  \"iterations\": {d},\n", .{measured_iterations});
    try writer.writeAll("  \"benchmarks\": [\n");

    for (results, 0..) |r, i| {
        try writer.print(
            "    {{\"name\": \"{s}\", \"median_ns\": {d}, \"min_ns\": {d}, \"max_ns\": {d}}}",
            .{ r.name, r.median_ns, r.min_ns, r.max_ns },
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

fn findBaseline(baselines: []const Baseline, name: []const u8) ?u64 {
    for (baselines) |b| {
        if (std.mem.eql(u8, b.name, name)) return b.median_ns;
    }
    return null;
}

fn checkRegressions(results: []const BenchmarkResult, baselines: []const Baseline) bool {
    var any_regression = false;
    const stderr = stderrWriter();

    for (results) |r| {
        const base = findBaseline(baselines, r.name) orelse continue;
        const ratio = @as(f64, @floatFromInt(r.median_ns)) / @as(f64, @floatFromInt(base));
        if (ratio > 1.0 + regression_threshold) {
            stderr.print("  REGRESSION: {s} — {d:.1}% slower (baseline: ", .{
                r.name, (ratio - 1.0) * 100.0,
            }) catch {};
            printDuration(stderr, base) catch {};
            stderr.writeAll(", current: ") catch {};
            printDuration(stderr, r.median_ns) catch {};
            stderr.writeAll(")\n") catch {};
            any_regression = true;
        }
    }

    return any_regression;
}

// ── Main ────────────────────────────────────────────────────────────────

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse CLI flags.
    var check_mode = false;
    var update_mode = false;
    var json_mode = false;

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);
    for (args[1..]) |arg| {
        if (std.mem.eql(u8, arg, "--check")) check_mode = true;
        if (std.mem.eql(u8, arg, "--update")) update_mode = true;
        if (std.mem.eql(u8, arg, "--json")) json_mode = true;
    }

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

    // Run all benchmarks.
    var results: [all_benchmarks.len]BenchmarkResult = undefined;
    for (all_benchmarks, 0..) |def, i| {
        try stderr.print("  Running {s}...\n", .{def.name});
        results[i] = runBenchmark(def, &ctx);
    }

    // Load baselines if they exist.
    var parsed_baselines = try readBaselines(allocator);
    defer if (parsed_baselines) |*p| p.deinit();

    const baseline_list: ?[]const Baseline = if (parsed_baselines) |p| p.value.benchmarks else null;

    // Output results.
    if (json_mode) {
        try printJson(stdoutWriter(), &results);
    }
    try printTable(stderr, &results, baseline_list);

    // Update baselines if requested.
    if (update_mode) {
        try writeBaselines(&results);
        try stderr.writeAll("  Baselines updated in bench/baselines.json\n\n");
    }

    // Check for regressions if requested.
    if (check_mode) {
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
}
