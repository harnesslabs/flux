//! Explicit FEEC operator context.
//!
//! This public family-specific context exposes FEEC/Whitney-based operators.

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

        pub fn hodgeStar(
            self: *Self,
            comptime k: comptime_int,
        ) !*const @import("hodge_star.zig").AssembledHodgeStar(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            return self.base.hodgeStar(k);
        }

        pub fn hodgeStarInverse(
            self: *Self,
            comptime k: comptime_int,
        ) !*const @import("hodge_star.zig").AssembledHodgeStarInverse(
            cochain.Cochain(MeshType, MeshType.topological_dimension - k, cochain.Dual),
        ) {
            return self.base.hodgeStarInverse(k);
        }

        pub fn codifferential(
            self: *Self,
            comptime k: comptime_int,
        ) !*const @import("codifferential.zig").AssembledCodifferential(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            return self.base.codifferential(k);
        }

        pub fn laplacian(
            self: *Self,
            comptime k: comptime_int,
        ) !*const @import("laplacian.zig").AssembledLaplacian(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            return self.base.laplacian(k);
        }
    };
}
