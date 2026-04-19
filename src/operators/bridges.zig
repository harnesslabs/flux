//! Explicit DEC/FEEC bridge operators.
//!
//! `WhitneyInterpolation` lifts simplicial cochain storage into a Whitney FEEC
//! space view, and `DeRhamProjection` returns Whitney-form coefficients back to
//! simplicial storage. For current lowest-order Whitney spaces, the same bridge
//! also provides the de Rham map from a continuous form into simplicial
//! coefficients by evaluating the simplex DOFs associated with vertices, edges,
//! and faces.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const feec = @import("../forms/feec.zig");
const topology = @import("../topology/mesh.zig");

pub fn WhitneyInterpolation(comptime MeshType: type, comptime k: comptime_int) type {
    const SpaceType = feec.WhitneySpace(MeshType, k);
    const Storage = cochain.Cochain(MeshType, k, cochain.Primal);

    return struct {
        const Self = @This();

        space: SpaceType,

        pub fn init(mesh: *const MeshType) Self {
            return .{ .space = SpaceType.init(mesh) };
        }

        pub fn apply(self: Self, coefficients: *const Storage) feec.Form(SpaceType) {
            return self.space.view(coefficients);
        }
    };
}

pub fn DeRhamProjection(comptime SpaceType: type) type {
    const Storage = SpaceType.StorageT;

    return struct {
        const Self = @This();

        space: SpaceType,

        pub fn init(space: SpaceType) Self {
            return .{ .space = space };
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, form: feec.Form(SpaceType)) !Storage {
            std.debug.assert(form.space.mesh == self.space.mesh);
            return form.intoCoefficients(allocator);
        }

        pub fn project(self: Self, allocator: std.mem.Allocator, continuous_form: anytype) !Storage {
            comptime {
                if (SpaceType.family != feec.Whitney) {
                    @compileError("DeRhamProjection.project currently supports only Whitney spaces");
                }
                if (Storage.duality != cochain.Primal) {
                    @compileError("DeRhamProjection.project currently supports only primal cochains");
                }
            }

            var coefficients = try Storage.init(allocator, self.space.mesh);
            errdefer coefficients.deinit(allocator);
            self.projectInto(&coefficients, continuous_form);
            return coefficients;
        }

        pub fn projectInto(self: Self, coefficients: *Storage, continuous_form: anytype) void {
            comptime {
                if (SpaceType.family != feec.Whitney) {
                    @compileError("DeRhamProjection.projectInto currently supports only Whitney spaces");
                }
                if (Storage.duality != cochain.Primal) {
                    @compileError("DeRhamProjection.projectInto currently supports only primal cochains");
                }
            }

            std.debug.assert(coefficients.mesh == self.space.mesh);
            projectByPrimalSimplexPairing(SpaceType.degree, self.space.mesh, coefficients.values, continuous_form);
        }
    };
}

fn projectByPrimalSimplexPairing(
    comptime simplex_degree: comptime_int,
    mesh: anytype,
    values: []f64,
    continuous_form: anytype,
) void {
    const coords = mesh.vertices.slice().items(.coords);
    switch (simplex_degree) {
        0 => {
            for (values, coords) |*value, point| {
                value.* = pairWithPrimalSimplex(0, continuous_form, .{point});
            }
        },
        1 => {
            const simplex_vertices = mesh.simplices(1).items(.vertices);
            for (values, simplex_vertices) |*value, simplex| {
                value.* = pairWithPrimalSimplex(1, continuous_form, .{
                    coords[simplex[0]],
                    coords[simplex[1]],
                });
            }
        },
        2 => {
            const simplex_vertices = mesh.simplices(2).items(.vertices);
            for (values, simplex_vertices) |*value, simplex| {
                value.* = pairWithPrimalSimplex(2, continuous_form, .{
                    coords[simplex[0]],
                    coords[simplex[1]],
                    coords[simplex[2]],
                });
            }
        },
        else => @compileError("DeRhamProjection.projectInto currently supports only degrees 0, 1, and 2"),
    }
}

fn pairWithPrimalSimplex(
    comptime simplex_degree: comptime_int,
    continuous_form: anytype,
    simplex_vertices: anytype,
) f64 {
    const ContinuousForm = @TypeOf(continuous_form);

    switch (simplex_degree) {
        0 => {
            const point = simplex_vertices[0];
            if (@hasDecl(ContinuousForm, "pointValue")) return continuous_form.pointValue(point);
            if (@hasDecl(ContinuousForm, "evaluate")) return continuous_form.evaluate(point);
            @compileError("continuous 0-form must declare `pointValue(point)` or `evaluate(point)`");
        },
        1 => {
            const p0 = simplex_vertices[0];
            const p1 = simplex_vertices[1];
            if (@hasDecl(ContinuousForm, "integrateEdge")) return continuous_form.integrateEdge(p0, p1);
            if (!@hasDecl(ContinuousForm, "evaluate")) {
                @compileError("continuous 1-form must declare `integrateEdge(p0, p1)` or `evaluate(point)`");
            }

            const nodes = [_]f64{
                0.5 * (1.0 - std.math.sqrt(3.0 / 5.0)),
                0.5,
                0.5 * (1.0 + std.math.sqrt(3.0 / 5.0)),
            };
            const weights = [_]f64{ 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 };

            var tangent: @TypeOf(p0) = undefined;
            inline for (0..p0.len) |axis| {
                tangent[axis] = p1[axis] - p0[axis];
            }

            var sum: f64 = 0.0;
            for (nodes, weights) |node, weight| {
                var point: @TypeOf(p0) = undefined;
                inline for (0..p0.len) |axis| {
                    point[axis] = (1.0 - node) * p0[axis] + node * p1[axis];
                }
                sum += weight * innerProduct(continuous_form.evaluate(point), tangent);
            }
            return sum;
        },
        2 => {
            const p0 = simplex_vertices[0];
            const p1 = simplex_vertices[1];
            const p2 = simplex_vertices[2];
            if (@hasDecl(ContinuousForm, "integrateFace")) return continuous_form.integrateFace(p0, p1, p2);
            if (!@hasDecl(ContinuousForm, "evaluate")) {
                @compileError("continuous 2-form must declare `integrateFace(p0, p1, p2)` or `evaluate(point)`");
            }

            const barycentric_nodes = [_][3]f64{
                .{ 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 },
                .{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 },
                .{ 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 },
            };

            if (p0.len == 2) {
                const signed_area_twice =
                    (p1[0] - p0[0]) * (p2[1] - p0[1]) -
                    (p1[1] - p0[1]) * (p2[0] - p0[0]);

                var sum: f64 = 0.0;
                for (barycentric_nodes) |lambda| {
                    const point = [2]f64{
                        lambda[0] * p0[0] + lambda[1] * p1[0] + lambda[2] * p2[0],
                        lambda[0] * p0[1] + lambda[1] * p1[1] + lambda[2] * p2[1],
                    };
                    sum += continuous_form.evaluate(point);
                }
                return (signed_area_twice / 6.0) * sum;
            }

            if (p0.len == 3) {
                const oriented_area_pseudovector_twice = triangleAreaPseudovectorTwice3(p0, p1, p2);

                var sum: f64 = 0.0;
                for (barycentric_nodes) |lambda| {
                    const point = [3]f64{
                        lambda[0] * p0[0] + lambda[1] * p1[0] + lambda[2] * p2[0],
                        lambda[0] * p0[1] + lambda[1] * p1[1] + lambda[2] * p2[1],
                        lambda[0] * p0[2] + lambda[1] * p1[2] + lambda[2] * p2[2],
                    };
                    sum += innerProduct(continuous_form.evaluate(point), oriented_area_pseudovector_twice);
                }
                return sum / 6.0;
            }

            @compileError("DeRhamProjection.projectInto degree-2 quadrature currently supports only 2D and 3D embeddings");
        },
        else => @compileError("DeRhamProjection.projectInto currently supports only degrees 0, 1, and 2"),
    }
}

fn triangleAreaPseudovectorTwice3(p0: [3]f64, p1: [3]f64, p2: [3]f64) [3]f64 {
    const edge_a = [3]f64{
        p1[0] - p0[0],
        p1[1] - p0[1],
        p1[2] - p0[2],
    };
    const edge_b = [3]f64{
        p2[0] - p0[0],
        p2[1] - p0[1],
        p2[2] - p0[2],
    };
    return .{
        edge_a[1] * edge_b[2] - edge_a[2] * edge_b[1],
        edge_a[2] * edge_b[0] - edge_a[0] * edge_b[2],
        edge_a[0] * edge_b[1] - edge_a[1] * edge_b[0],
    };
}

fn innerProduct(lhs: anytype, rhs: @TypeOf(lhs)) f64 {
    var sum: f64 = 0.0;
    inline for (0..lhs.len) |axis| {
        sum += lhs[axis] * rhs[axis];
    }
    return sum;
}

const Mesh2D = topology.Mesh(2, 2);
const C1 = cochain.Cochain(Mesh2D, 1, cochain.Primal);

test "Whitney interpolation borrows cochain storage into FEEC space" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var coefficients = try C1.init(allocator, &mesh);
    defer coefficients.deinit(allocator);
    coefficients.values[0] = 3.5;

    const interpolate = WhitneyInterpolation(Mesh2D, 1).init(&mesh);
    const form = interpolate.apply(&coefficients);

    try testing.expect(form.coefficientsConst().values.ptr == coefficients.values.ptr);
    try testing.expectEqual(@as(f64, 3.5), form.coefficientsConst().values[0]);
}

test "de Rham projection preserves coefficients for Whitney forms" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var coefficients = try C1.init(allocator, &mesh);
    defer coefficients.deinit(allocator);
    for (coefficients.values, 0..) |*value, i| {
        value.* = @floatFromInt(i + 1);
    }

    const Space1 = feec.WhitneySpace(Mesh2D, 1);
    const interpolate = WhitneyInterpolation(Mesh2D, 1).init(&mesh);
    const project = DeRhamProjection(Space1).init(Space1.init(&mesh));

    const form = interpolate.apply(&coefficients);
    var projected = try project.apply(allocator, form);
    defer projected.deinit(allocator);

    try testing.expectEqualSlices(f64, coefficients.values, projected.values);
}

test "de Rham projection maps a continuous constant one-form onto primal 1-cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const Space1 = feec.WhitneySpace(Mesh2D, 1);
    const project = DeRhamProjection(Space1).init(Space1.init(&mesh));

    const ConstantOneForm = struct {
        pub fn evaluate(_: @This(), point: [2]f64) [2]f64 {
            _ = point;
            return .{ 2.0, -1.0 };
        }
    };

    var projected = try project.project(allocator, ConstantOneForm{});
    defer projected.deinit(allocator);

    const edge_vertices = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (projected.values, edge_vertices) |value, edge| {
        const p0 = coords[edge[0]];
        const p1 = coords[edge[1]];
        const tangent = [2]f64{ p1[0] - p0[0], p1[1] - p0[1] };
        try testing.expectApproxEqAbs(innerProduct([2]f64{ 2.0, -1.0 }, tangent), value, 1e-12);
    }
}

test "de Rham projection uses continuous edge DOFs when provided" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const Space1 = feec.WhitneySpace(Mesh2D, 1);
    const project = DeRhamProjection(Space1).init(Space1.init(&mesh));

    const ExactQuadraticOneForm = struct {
        pub fn evaluate(_: @This(), point: [2]f64) [2]f64 {
            return .{ point[0] * point[0], 0.0 };
        }

        pub fn integrateEdge(_: @This(), p0: [2]f64, p1: [2]f64) f64 {
            const dx = p1[0] - p0[0];
            return dx * (p0[0] * p0[0] + p0[0] * dx + (dx * dx) / 3.0);
        }
    };

    var projected = try project.project(allocator, ExactQuadraticOneForm{});
    defer projected.deinit(allocator);

    const edge_vertices = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (projected.values, edge_vertices) |value, edge| {
        const p0 = coords[edge[0]];
        const p1 = coords[edge[1]];
        const dx = p1[0] - p0[0];
        const exact = dx * (p0[0] * p0[0] + p0[0] * dx + (dx * dx) / 3.0);
        try testing.expectApproxEqAbs(exact, value, 1e-12);
    }
}

test "de Rham projection maps a continuous constant two-form onto primal 2-cochains" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const Space2 = feec.WhitneySpace(Mesh2D, 2);
    const project = DeRhamProjection(Space2).init(Space2.init(&mesh));

    const ConstantTwoForm = struct {
        pub fn evaluate(_: @This(), point: [2]f64) f64 {
            _ = point;
            return 2.0;
        }
    };

    var projected = try project.project(allocator, ConstantTwoForm{});
    defer projected.deinit(allocator);

    for (projected.values) |value| {
        try testing.expectApproxEqAbs(1.0, value, 1e-12);
    }
}
