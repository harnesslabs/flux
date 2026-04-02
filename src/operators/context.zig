//! Operator assembly context bound to a single mesh.
//!
//! The mesh owns topology and geometric primitives. The operator context owns
//! only the assembled operator state actually requested by a given problem.

const std = @import("std");
const cochain = @import("../forms/cochain.zig");
const laplacian_mod = @import("laplacian.zig");

fn LaplacianSlotsType(comptime MeshType: type) type {
    const count = MeshType.topological_dimension + 1;
    var fields: [count]std.builtin.Type.StructField = undefined;

    for (0..count) |i| {
        const InputType = cochain.Cochain(MeshType, i, cochain.Primal);
        const SlotType = ?laplacian_mod.AssembledLaplacian(InputType);
        fields[i] = .{
            .name = std.fmt.comptimePrint("{d}", .{i}),
            .type = SlotType,
            .default_value_ptr = null,
            .is_comptime = false,
            .alignment = @alignOf(SlotType),
        };
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &fields,
        .decls = &.{},
        .is_tuple = true,
    } });
}

/// Per-mesh owner of assembled DEC operators.
///
/// The context is the explicit efficient path: callers request only the
/// operators they need, then reuse them for repeated application.
pub fn OperatorContext(comptime MeshType: type) type {
    const LaplacianSlots = LaplacianSlotsType(MeshType);

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        mesh: *const MeshType,
        laplacians: LaplacianSlots,

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) Self {
            return .{
                .allocator = allocator,
                .mesh = mesh,
                .laplacians = std.mem.zeroes(LaplacianSlots),
            };
        }

        pub fn deinit(self: *Self) void {
            inline for (0..MeshType.topological_dimension + 1) |k| {
                if (self.laplacians[k]) |*op| {
                    op.deinit(self.allocator);
                }
            }
        }

        /// Ensure Δₖ is assembled and available for repeated application.
        pub fn withLaplacian(self: *Self, comptime k: comptime_int) !void {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no Laplacian degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            if (self.laplacians[k] == null) {
                self.laplacians[k] = try laplacian_mod.assemble_for_degree(MeshType, k, self.allocator, self.mesh);
            }
        }

        /// Access an assembled Δₖ. The operator must already have been requested
        /// via `withLaplacian(k)`.
        pub fn laplacian(self: *Self, comptime k: comptime_int) *const laplacian_mod.AssembledLaplacian(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no Laplacian degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }
            std.debug.assert(self.laplacians[k] != null);
            return &self.laplacians[k].?;
        }
    };
}
