//! FEEC form-space abstractions layered over shared cochain storage.

const std = @import("std");
const cochain = @import("cochain.zig");

pub const Whitney = struct {};

pub fn Space(comptime MeshType: type, comptime k: comptime_int, comptime Family: type) type {
    const Storage = cochain.Cochain(MeshType, k, cochain.Primal);

    return struct {
        const Self = @This();

        pub const MeshT = MeshType;
        pub const degree = k;
        pub const family = Family;
        pub const StorageT = Storage;

        mesh: *const MeshType,

        pub fn init(mesh: *const MeshType) Self {
            return .{ .mesh = mesh };
        }

        pub fn view(self: Self, coefficients: *const Storage) Form(Self) {
            std.debug.assert(coefficients.mesh == self.mesh);
            return .{
                .space = self,
                .storage = .{ .borrowed = coefficients },
            };
        }
    };
}

pub fn WhitneySpace(comptime MeshType: type, comptime k: comptime_int) type {
    return Space(MeshType, k, Whitney);
}

pub fn Form(comptime SpaceType: type) type {
    const Storage = SpaceType.StorageT;

    return struct {
        const Self = @This();

        pub const SpaceT = SpaceType;
        pub const MeshT = SpaceType.MeshT;
        pub const degree = SpaceType.degree;

        pub const StorageHandle = union(enum) {
            borrowed: *const Storage,
            owned: Storage,
        };

        space: SpaceType,
        storage: StorageHandle,

        pub fn initOwned(allocator: std.mem.Allocator, space: SpaceType) !Self {
            return .{
                .space = space,
                .storage = .{ .owned = try Storage.init(allocator, space.mesh) },
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            switch (self.storage) {
                .borrowed => {},
                .owned => |*owned| owned.deinit(allocator),
            }
        }

        pub fn coefficientsConst(self: *const Self) *const Storage {
            return switch (self.storage) {
                .borrowed => |borrowed| borrowed,
                .owned => |*owned| owned,
            };
        }

        pub fn coefficientsMut(self: *Self) *Storage {
            return switch (self.storage) {
                .borrowed => @panic("cannot mutably borrow coefficients from a FEEC view"),
                .owned => |*owned| owned,
            };
        }
    };
}
