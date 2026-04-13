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
    };
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
