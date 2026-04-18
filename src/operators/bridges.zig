//! Explicit DEC/FEEC bridge operators.
//!
//! `WhitneyInterpolation` lifts simplicial cochain storage into a Whitney FEEC
//! space view, and `DeRhamProjection` returns FEEC form coefficients back to
//! simplicial storage.

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
            std.debug.assert(coefficients.mesh == self.space.mesh);
            switch (SpaceType.degree) {
                0 => projectWhitneyZero(self.space.mesh, coefficients.values, continuous_form),
                1 => projectWhitneyOne(self.space.mesh, coefficients.values, continuous_form),
                2 => projectWhitneyTwo(self.space.mesh, coefficients.values, continuous_form),
                else => @compileError("DeRhamProjection.project currently supports only degrees 0, 1, and 2"),
            }
        }
    };
}

fn projectWhitneyZero(mesh: anytype, values: []f64, continuous_form: anytype) void {
    const coords = mesh.vertices.slice().items(.coords);
    for (values, coords) |*value, point| {
        value.* = scalarValue(continuous_form.evaluate(point));
    }
}

fn projectWhitneyOne(mesh: anytype, values: []f64, continuous_form: anytype) void {
    const edge_vertices = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_vertices) |*value, edge| {
        const p0 = coords[edge[0]];
        const p1 = coords[edge[1]];
        var midpoint: @TypeOf(p0) = undefined;
        var tangent: @TypeOf(p0) = undefined;
        inline for (0..p0.len) |axis| {
            midpoint[axis] = 0.5 * (p0[axis] + p1[axis]);
            tangent[axis] = p1[axis] - p0[axis];
        }
        value.* = dot(continuous_form.evaluate(midpoint), tangent);
    }
}

fn projectWhitneyTwo(mesh: anytype, values: []f64, continuous_form: anytype) void {
    const face_vertices = mesh.simplices(2).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    const face_areas = mesh.simplices(2).items(.volume);

    switch (@TypeOf(mesh.*).topological_dimension) {
        2 => {
            for (values, face_vertices, face_areas) |*value, face, area| {
                const centroid = triangleCentroid(coords[face[0]], coords[face[1]], coords[face[2]]);
                value.* = scalarValue(continuous_form.evaluate(centroid)) * area;
            }
        },
        3 => {
            for (values, face_vertices) |*value, face| {
                const p0 = coords[face[0]];
                const p1 = coords[face[1]];
                const p2 = coords[face[2]];
                const centroid = triangleCentroid(p0, p1, p2);
                value.* = dot(continuous_form.evaluate(centroid), triangleAreaVector(p0, p1, p2));
            }
        },
        else => @compileError("DeRhamProjection.project degree-2 currently supports only 2D and 3D meshes"),
    }
}

fn triangleCentroid(p0: anytype, p1: @TypeOf(p0), p2: @TypeOf(p0)) @TypeOf(p0) {
    var centroid: @TypeOf(p0) = undefined;
    inline for (0..p0.len) |axis| {
        centroid[axis] = (p0[axis] + p1[axis] + p2[axis]) / 3.0;
    }
    return centroid;
}

fn triangleAreaVector(p0: [3]f64, p1: [3]f64, p2: [3]f64) [3]f64 {
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
        0.5 * (edge_a[1] * edge_b[2] - edge_a[2] * edge_b[1]),
        0.5 * (edge_a[2] * edge_b[0] - edge_a[0] * edge_b[2]),
        0.5 * (edge_a[0] * edge_b[1] - edge_a[1] * edge_b[0]),
    };
}

fn dot(lhs: anytype, rhs: @TypeOf(lhs)) f64 {
    var sum: f64 = 0.0;
    inline for (0..lhs.len) |axis| {
        sum += lhs[axis] * rhs[axis];
    }
    return sum;
}

fn scalarValue(value: anytype) f64 {
    return value;
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

test "de Rham projection samples a constant one-form onto primal 1-cochains" {
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
        try testing.expectEqual(dot([2]f64{ 2.0, -1.0 }, tangent), value);
    }
}

test "de Rham projection samples a constant two-form onto primal 2-cochains in 2D" {
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
        try testing.expectEqual(@as(f64, 1.0), value);
    }
}
