//! Explicit DEC operator context.
//!
//! This public family-specific context exposes only DEC/topological operators.

const std = @import("std");
const base_context = @import("context.zig");
const cochain = @import("../forms/cochain.zig");

pub fn OperatorContext(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        mesh: *const MeshType,
        base: *base_context.OperatorContext(MeshType),

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !*Self {
            const self = try allocator.create(Self);
            errdefer allocator.destroy(self);

            self.* = .{
                .allocator = allocator,
                .mesh = mesh,
                .base = try base_context.OperatorContext(MeshType).init(allocator, mesh),
            };
            return self;
        }

        pub fn deinit(self: *Self) void {
            self.base.deinit();
            self.allocator.destroy(self);
        }

        pub fn exteriorDerivative(
            self: *Self,
            comptime Duality: type,
            comptime k: comptime_int,
        ) !*const @import("exterior_derivative.zig").AssembledExteriorDerivative(
            cochain.Cochain(MeshType, k, Duality),
        ) {
            return self.base.exteriorDerivative(Duality, k);
        }
    };
}
