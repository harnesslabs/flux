//! Discrete wedge product on primal simplicial cochains.
//!
//! The current implementation targets the lowest-order FEEC/Whitney setting:
//! cochains are interpreted as Whitney form coefficients, the smooth wedge is
//! taken in that space, and the result is projected back to a cochain by the
//! de Rham map. We compute the induced local simplex formula directly rather
//! than materializing the interpolated forms.
//!
//! Algebraic properties:
//!   - graded commutativity holds exactly;
//!   - Leibniz with the exterior derivative holds exactly;
//!   - associativity holds for closed forms, but not in general.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const feec = @import("../forms/feec.zig");
const topology = @import("../topology/mesh.zig");
const bridges = @import("bridges.zig");
const exterior_derivative = @import("exterior_derivative.zig");

fn WedgeResult(comptime LeftType: type, comptime RightType: type) type {
    return feec.Form(feec.WhitneySpace(LeftType.MeshT, LeftType.degree + RightType.degree));
}

fn validateWedgeInputs(comptime LeftType: type, comptime RightType: type) void {
    if (!@hasDecl(LeftType, "SpaceT") or !@hasDecl(RightType, "SpaceT")) {
        @compileError("wedge requires FEEC form inputs rather than bare cochains");
    }
    if (LeftType.MeshT != RightType.MeshT) {
        @compileError("wedge requires both FEEC forms to use the same mesh type");
    }
    if (LeftType.SpaceT.family != feec.Whitney or RightType.SpaceT.family != feec.Whitney) {
        @compileError("wedge currently supports only Whitney FEEC spaces");
    }
    if (LeftType.degree + RightType.degree > LeftType.MeshT.topological_dimension) {
        @compileError(std.fmt.comptimePrint(
            "wedge degree {d} + {d} exceeds mesh topological dimension {d}",
            .{ LeftType.degree, RightType.degree, LeftType.MeshT.topological_dimension },
        ));
    }
}

pub fn wedge(
    allocator: std.mem.Allocator,
    left: anytype,
    right: anytype,
) !WedgeResult(@TypeOf(left), @TypeOf(right)) {
    const LeftType = @TypeOf(left);
    const RightType = @TypeOf(right);
    comptime validateWedgeInputs(LeftType, RightType);

    const OutputType = WedgeResult(LeftType, RightType);
    std.debug.assert(left.space.mesh == right.space.mesh);

    var output = try OutputType.initOwned(
        allocator,
        feec.WhitneySpace(LeftType.MeshT, LeftType.degree + RightType.degree).init(left.space.mesh),
    );
    errdefer output.deinit(allocator);

    try wedgeCoefficientsInto(
        output.coefficientsMut().values,
        left.coefficientsConst().*,
        right.coefficientsConst().*,
    );
    return output;
}

fn wedgeCoefficientsInto(output_values: []f64, left: anytype, right: anytype) !void {
    const LeftType = @TypeOf(left);
    const RightType = @TypeOf(right);
    std.debug.assert(left.mesh == right.mesh);

    const degree_left = LeftType.degree;
    const degree_right = RightType.degree;
    const degree_out = degree_left + degree_right;

    if (degree_out == 0) {
        for (output_values, left.values, right.values) |*out, left_value, right_value| {
            out.* = left_value * right_value;
        }
        return;
    }

    const left_values = left.values;
    const right_values = right.values;
    const output_simplices = left.mesh.simplices(degree_out).items(.vertices);

    const coefficient = wedgeCoefficient(degree_left, degree_right);

    for (output_simplices, 0..) |simplex_vertices, out_index| {
        var sum: f64 = 0.0;
        const permutation_count = factorial(degree_out + 1);
        var permutation_index: usize = 0;
        while (permutation_index < permutation_count) : (permutation_index += 1) {
            const permutation = nthPermutation(degree_out + 1, permutation_index);

            var left_vertices: [degree_left + 1]u32 = undefined;
            if (degree_left == 0) {
                left_vertices[0] = simplex_vertices[permutation[0]];
            } else {
                for (0..degree_left + 1) |i| {
                    left_vertices[i] = simplex_vertices[permutation[i]];
                }
            }

            var right_vertices: [degree_right + 1]u32 = undefined;
            if (degree_right == 0) {
                right_vertices[0] = simplex_vertices[permutation[degree_left]];
            } else {
                for (0..degree_right + 1) |i| {
                    right_vertices[i] = simplex_vertices[permutation[degree_left + i]];
                }
            }

            const left_value = if (degree_left == 0)
                left.values[left_vertices[0]]
            else blk: {
                const left_simplices = left.mesh.simplices(degree_left).items(.vertices);
                const simplex_index = findSimplexIndex(left_vertices.len, left_simplices, left_vertices) orelse continue;
                const orientation = orientationSign(left_vertices.len, left_vertices, left_simplices[simplex_index]);
                break :blk @as(f64, @floatFromInt(orientation)) * left_values[simplex_index];
            };

            const right_value = if (degree_right == 0)
                right.values[right_vertices[0]]
            else blk: {
                const right_simplices = right.mesh.simplices(degree_right).items(.vertices);
                const simplex_index = findSimplexIndex(right_vertices.len, right_simplices, right_vertices) orelse continue;
                const orientation = orientationSign(right_vertices.len, right_vertices, right_simplices[simplex_index]);
                break :blk @as(f64, @floatFromInt(orientation)) * right_values[simplex_index];
            };

            sum += @as(f64, @floatFromInt(permutationSign(degree_out + 1, permutation))) * left_value * right_value;
        }
        output_values[out_index] = coefficient * sum;
    }
}

fn wedgeCoefficient(comptime degree_left: comptime_int, comptime degree_right: comptime_int) f64 {
    return 1.0 / @as(f64, @floatFromInt(factorial(degree_left + degree_right + 1)));
}

fn factorial(comptime n: comptime_int) comptime_int {
    var result: comptime_int = 1;
    inline for (1..n + 1) |i| {
        result *= i;
    }
    return result;
}

fn findSimplexIndex(
    comptime len: comptime_int,
    simplices: []const [len]u32,
    vertices: [len]u32,
) ?usize {
    const key = canonicalizeVertices(len, vertices);
    for (simplices, 0..) |simplex_vertices, index| {
        if (std.mem.eql(u32, &canonicalizeVertices(len, simplex_vertices), &key)) {
            return index;
        }
    }
    return null;
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
    for (oriented_vertices, 0..) |vertex, oriented_index| {
        var found = false;
        for (global_vertices, 0..) |global_vertex, global_index| {
            if (used[global_index]) continue;
            if (vertex != global_vertex) continue;
            permutation[oriented_index] = @intCast(global_index);
            used[global_index] = true;
            found = true;
            break;
        }
        std.debug.assert(found);
    }
    return permutationSign(len, permutation);
}

fn permutationSign(comptime len: comptime_int, permutation: [len]u8) i8 {
    var inversions: usize = 0;
    inline for (0..len) |i| {
        inline for (i + 1..len) |j| {
            if (permutation[i] > permutation[j]) {
                inversions += 1;
            }
        }
    }
    return if (inversions % 2 == 0) 1 else -1;
}

fn nthPermutation(comptime len: comptime_int, index: usize) [len]u8 {
    var available: [len]u8 = undefined;
    for (0..len) |i| {
        available[i] = @intCast(i);
    }

    var result: [len]u8 = undefined;
    var remaining = index;
    var available_len: usize = len;
    inline for (0..len) |position| {
        const block_size = factorial(len - 1 - position);
        const choice_index = remaining / block_size;
        remaining %= block_size;

        result[position] = available[choice_index];
        var i = choice_index;
        while (i + 1 < available_len) : (i += 1) {
            available[i] = available[i + 1];
        }
        available_len -= 1;
    }
    return result;
}

const Mesh2D = topology.Mesh(2, 2);
const Mesh3D = topology.Mesh(3, 3);
const C0 = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const C1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const C2 = cochain.Cochain(Mesh2D, 2, cochain.Primal);
const C3D0 = cochain.Cochain(Mesh3D, 0, cochain.Primal);
const C3D1 = cochain.Cochain(Mesh3D, 1, cochain.Primal);
const C3D2 = cochain.Cochain(Mesh3D, 2, cochain.Primal);
const C3D3 = cochain.Cochain(Mesh3D, 3, cochain.Primal);

test "compile-time: wedge degree arithmetic yields the sum degree" {
    comptime {
        const Whitney2D1 = feec.Form(feec.WhitneySpace(Mesh2D, 1));
        const Result2D = WedgeResult(Whitney2D1, Whitney2D1);
        if (Result2D.SpaceT.degree != 2) {
            @compileError("2D FEEC wedge should produce a degree-2 FEEC form");
        }

        const Whitney3D1 = feec.Form(feec.WhitneySpace(Mesh3D, 1));
        const Whitney3D2 = feec.Form(feec.WhitneySpace(Mesh3D, 2));
        const Result3D = WedgeResult(Whitney3D1, Whitney3D2);
        if (Result3D.SpaceT.degree != 3) {
            @compileError("3D FEEC wedge should produce a degree-3 FEEC form");
        }
    }
}

test "wedge of 0-forms is pointwise multiplication" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var alpha = try C0.init(allocator, &mesh);
    defer alpha.deinit(allocator);
    var beta = try C0.init(allocator, &mesh);
    defer beta.deinit(allocator);

    for (alpha.values, 0..) |*value, i| {
        value.* = @floatFromInt(i + 1);
    }
    for (beta.values, 0..) |*value, i| {
        value.* = 2.0 * @as(f64, @floatFromInt(i + 3));
    }

    const interpolate = bridges.WhitneyInterpolation(Mesh2D, 0).init(&mesh);
    var product = try wedge(allocator, interpolate.apply(&alpha), interpolate.apply(&beta));
    defer product.deinit(allocator);

    for (product.coefficientsConst().values, 0..) |value, i| {
        try testing.expectApproxEqAbs(alpha.values[i] * beta.values[i], value, 1e-15);
    }
}

test "wedge operates on FEEC Whitney forms instead of bare cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var alpha_coefficients = try C1.init(allocator, &mesh);
    defer alpha_coefficients.deinit(allocator);
    var beta_coefficients = try C1.init(allocator, &mesh);
    defer beta_coefficients.deinit(allocator);

    for (alpha_coefficients.values, 0..) |*value, i| {
        value.* = @floatFromInt(i + 1);
    }
    for (beta_coefficients.values, 0..) |*value, i| {
        value.* = @floatFromInt(2 * i + 1);
    }

    const interpolate = bridges.WhitneyInterpolation(Mesh2D, 1).init(&mesh);
    const alpha = interpolate.apply(&alpha_coefficients);
    const beta = interpolate.apply(&beta_coefficients);

    var product = try wedge(allocator, alpha, beta);
    defer product.deinit(allocator);

    try testing.expect(@TypeOf(product).SpaceT.degree == 2);
    try testing.expect(product.coefficientsConst().mesh == &mesh);
}

test "wedge of 0-form and 1-form averages the endpoint values" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var scalar = try C0.init(allocator, &mesh);
    defer scalar.deinit(allocator);
    const coords = mesh.vertices.slice().items(.coords);
    for (scalar.values, coords) |*value, coord| {
        value.* = coord[0] + 2.0 * coord[1];
    }

    var one_form = try C1.init(allocator, &mesh);
    defer one_form.deinit(allocator);
    for (one_form.values, 0..) |*value, i| {
        value.* = @as(f64, @floatFromInt(i + 1));
    }

    const interpolate_scalar = bridges.WhitneyInterpolation(Mesh2D, 0).init(&mesh);
    const interpolate_one_form = bridges.WhitneyInterpolation(Mesh2D, 1).init(&mesh);
    var product = try wedge(
        allocator,
        interpolate_scalar.apply(&scalar),
        interpolate_one_form.apply(&one_form),
    );
    defer product.deinit(allocator);

    const edges = mesh.simplices(1).items(.vertices);
    for (product.coefficientsConst().values, edges, one_form.values) |value, edge, edge_value| {
        const expected = 0.5 * (scalar.values[edge[0]] + scalar.values[edge[1]]) * edge_value;
        try testing.expectApproxEqAbs(expected, value, 1e-14);
    }
}

test "graded commutativity holds on random 3D 1- and 2-forms" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xA11CE_84_01);
    for (0..200) |_| {
        var alpha = try C3D1.init(allocator, &mesh);
        defer alpha.deinit(allocator);
        var beta = try C3D2.init(allocator, &mesh);
        defer beta.deinit(allocator);

        for (alpha.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;
        for (beta.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;

        const interpolate_one = bridges.WhitneyInterpolation(Mesh3D, 1).init(&mesh);
        const interpolate_two = bridges.WhitneyInterpolation(Mesh3D, 2).init(&mesh);
        var left = try wedge(allocator, interpolate_one.apply(&alpha), interpolate_two.apply(&beta));
        defer left.deinit(allocator);
        var right = try wedge(allocator, interpolate_two.apply(&beta), interpolate_one.apply(&alpha));
        defer right.deinit(allocator);

        for (left.coefficientsConst().values, right.coefficientsConst().values) |lhs, rhs| {
            try testing.expectApproxEqAbs(lhs, rhs, 1e-12);
        }
    }
}

test "associativity holds for random closed 1-forms on tetrahedral meshes" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 2, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xA11CE_84_02);
    for (0..100) |_| {
        var potential_a = try C3D0.init(allocator, &mesh);
        defer potential_a.deinit(allocator);
        var potential_b = try C3D0.init(allocator, &mesh);
        defer potential_b.deinit(allocator);
        var potential_c = try C3D0.init(allocator, &mesh);
        defer potential_c.deinit(allocator);

        for (potential_a.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;
        for (potential_b.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;
        for (potential_c.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;

        var alpha = try exterior_derivative.exterior_derivative(allocator, potential_a);
        defer alpha.deinit(allocator);
        var beta = try exterior_derivative.exterior_derivative(allocator, potential_b);
        defer beta.deinit(allocator);
        var gamma = try exterior_derivative.exterior_derivative(allocator, potential_c);
        defer gamma.deinit(allocator);

        const interpolate_one = bridges.WhitneyInterpolation(Mesh3D, 1).init(&mesh);
        var alpha_beta = try wedge(allocator, interpolate_one.apply(&alpha), interpolate_one.apply(&beta));
        defer alpha_beta.deinit(allocator);
        var left = try wedge(allocator, alpha_beta, interpolate_one.apply(&gamma));
        defer left.deinit(allocator);

        var beta_gamma = try wedge(allocator, interpolate_one.apply(&beta), interpolate_one.apply(&gamma));
        defer beta_gamma.deinit(allocator);
        var right = try wedge(allocator, interpolate_one.apply(&alpha), beta_gamma);
        defer right.deinit(allocator);

        for (left.coefficientsConst().values, right.coefficientsConst().values) |lhs, rhs| {
            try testing.expectApproxEqAbs(lhs, rhs, 1e-11);
        }
    }
}

test "Leibniz rule holds on random 2D 0- and 1-forms" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xA11CE_84_03);
    for (0..200) |_| {
        var alpha = try C0.init(allocator, &mesh);
        defer alpha.deinit(allocator);
        var beta = try C1.init(allocator, &mesh);
        defer beta.deinit(allocator);

        for (alpha.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;
        for (beta.values) |*value| value.* = rng.random().float(f64) * 2.0 - 1.0;

        const interpolate_zero = bridges.WhitneyInterpolation(Mesh2D, 0).init(&mesh);
        const interpolate_one = bridges.WhitneyInterpolation(Mesh2D, 1).init(&mesh);
        var alpha_wedge_beta = try wedge(allocator, interpolate_zero.apply(&alpha), interpolate_one.apply(&beta));
        defer alpha_wedge_beta.deinit(allocator);
        var left = try exterior_derivative.exterior_derivative(allocator, alpha_wedge_beta.coefficientsConst().*);
        defer left.deinit(allocator);

        var d_alpha = try exterior_derivative.exterior_derivative(allocator, alpha);
        defer d_alpha.deinit(allocator);
        var d_alpha_wedge_beta = try wedge(allocator, interpolate_one.apply(&d_alpha), interpolate_one.apply(&beta));
        defer d_alpha_wedge_beta.deinit(allocator);

        var d_beta = try exterior_derivative.exterior_derivative(allocator, beta);
        defer d_beta.deinit(allocator);
        var alpha_wedge_d_beta = try wedge(allocator, interpolate_zero.apply(&alpha), bridges.WhitneyInterpolation(Mesh2D, 2).init(&mesh).apply(&d_beta));
        defer alpha_wedge_d_beta.deinit(allocator);

        for (left.values, d_alpha_wedge_beta.coefficientsConst().values, alpha_wedge_d_beta.coefficientsConst().values) |lhs, rhs_a, rhs_b| {
            try testing.expectApproxEqAbs(lhs, rhs_a + rhs_b, 1e-11);
        }
    }
}
