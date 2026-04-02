//! Operator assembly context bound to a single mesh.
//!
//! The mesh owns topology and geometric primitives. The operator context owns
//! only the assembled operator state actually requested by a given problem.

const std = @import("std");
const cochain = @import("../forms/cochain.zig");
const codifferential_mod = @import("codifferential.zig");
const exterior_derivative_mod = @import("exterior_derivative.zig");
const hodge_star_mod = @import("hodge_star.zig");
const laplacian_mod = @import("laplacian.zig");

fn SlotsTupleType(
    comptime count: usize,
    comptime SlotTypeForIndex: fn (comptime usize) type,
) type {
    var fields: [count]std.builtin.Type.StructField = undefined;

    for (0..count) |i| {
        const SlotType = SlotTypeForIndex(i);
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

fn LaplacianSlotsType(comptime MeshType: type) type {
    const count = MeshType.topological_dimension + 1;
    return SlotsTupleType(count, struct {
        fn get(comptime i: usize) type {
            const InputType = cochain.Cochain(MeshType, i, cochain.Primal);
            return ?laplacian_mod.AssembledLaplacian(InputType);
        }
    }.get);
}

fn ExteriorDerivativeSlotsType(comptime MeshType: type, comptime Duality: type) type {
    const count = MeshType.topological_dimension;
    return SlotsTupleType(count, struct {
        fn get(comptime i: usize) type {
            const InputType = cochain.Cochain(MeshType, i, Duality);
            return ?exterior_derivative_mod.AssembledExteriorDerivative(InputType);
        }
    }.get);
}

fn HodgeStarSlotsType(comptime MeshType: type) type {
    const count = MeshType.topological_dimension + 1;
    return SlotsTupleType(count, struct {
        fn get(comptime i: usize) type {
            const InputType = cochain.Cochain(MeshType, i, cochain.Primal);
            return ?hodge_star_mod.AssembledHodgeStar(InputType);
        }
    }.get);
}

fn HodgeStarInverseSlotsType(comptime MeshType: type) type {
    const count = MeshType.topological_dimension + 1;
    return SlotsTupleType(count, struct {
        fn get(comptime i: usize) type {
            const InputType = cochain.Cochain(MeshType, MeshType.topological_dimension - i, cochain.Dual);
            return ?hodge_star_mod.AssembledHodgeStarInverse(InputType);
        }
    }.get);
}

fn CodifferentialSlotsType(comptime MeshType: type) type {
    const count = MeshType.topological_dimension;
    return SlotsTupleType(count, struct {
        fn get(comptime i: usize) type {
            const InputType = cochain.Cochain(MeshType, i + 1, cochain.Primal);
            return ?codifferential_mod.AssembledCodifferential(InputType);
        }
    }.get);
}

/// Per-mesh owner of assembled DEC operators.
///
/// The context is the explicit efficient path: callers request only the
/// operators they need, then reuse them for repeated application.
pub fn OperatorContext(comptime MeshType: type) type {
    const LaplacianSlots = LaplacianSlotsType(MeshType);
    const ExteriorDerivativePrimalSlots = ExteriorDerivativeSlotsType(MeshType, cochain.Primal);
    const ExteriorDerivativeDualSlots = ExteriorDerivativeSlotsType(MeshType, cochain.Dual);
    const HodgeStarSlots = HodgeStarSlotsType(MeshType);
    const HodgeStarInverseSlots = HodgeStarInverseSlotsType(MeshType);
    const CodifferentialSlots = CodifferentialSlotsType(MeshType);

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        mesh: *const MeshType,
        exterior_derivative_primal: ExteriorDerivativePrimalSlots,
        exterior_derivative_dual: ExteriorDerivativeDualSlots,
        hodge_stars: HodgeStarSlots,
        hodge_star_inverses: HodgeStarInverseSlots,
        codifferentials: CodifferentialSlots,
        laplacians: LaplacianSlots,

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) Self {
            return .{
                .allocator = allocator,
                .mesh = mesh,
                .exterior_derivative_primal = std.mem.zeroes(ExteriorDerivativePrimalSlots),
                .exterior_derivative_dual = std.mem.zeroes(ExteriorDerivativeDualSlots),
                .hodge_stars = std.mem.zeroes(HodgeStarSlots),
                .hodge_star_inverses = std.mem.zeroes(HodgeStarInverseSlots),
                .codifferentials = std.mem.zeroes(CodifferentialSlots),
                .laplacians = std.mem.zeroes(LaplacianSlots),
            };
        }

        pub fn deinit(self: *Self) void {
            inline for (0..MeshType.topological_dimension) |k| {
                if (self.exterior_derivative_primal[k]) |*op| op.deinit(self.allocator);
                if (self.exterior_derivative_dual[k]) |*op| op.deinit(self.allocator);
            }
            inline for (0..MeshType.topological_dimension + 1) |k| {
                if (self.hodge_stars[k]) |*op| op.deinit(self.allocator);
                if (self.hodge_star_inverses[k]) |*op| op.deinit(self.allocator);
            }
            inline for (0..MeshType.topological_dimension) |k| {
                if (self.codifferentials[k]) |*op| op.deinit(self.allocator);
            }
            inline for (0..MeshType.topological_dimension + 1) |k| {
                if (self.laplacians[k]) |*op| {
                    op.deinit(self.allocator);
                }
            }
        }

        pub fn withExteriorDerivative(self: *Self, comptime Duality: type, comptime k: comptime_int) !void {
            if (k < 0 or k >= MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no exterior derivative degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            if (Duality == cochain.Primal) {
                if (self.exterior_derivative_primal[k] == null) {
                    self.exterior_derivative_primal[k] = try exterior_derivative_mod.assemble_for_degree(
                        MeshType,
                        Duality,
                        k,
                        self.allocator,
                        self.mesh,
                    );
                }
                return;
            }

            if (Duality == cochain.Dual) {
                if (self.exterior_derivative_dual[k] == null) {
                    self.exterior_derivative_dual[k] = try exterior_derivative_mod.assemble_for_degree(
                        MeshType,
                        Duality,
                        k,
                        self.allocator,
                        self.mesh,
                    );
                }
                return;
            }

            @compileError("ExteriorDerivative duality must be cochain.Primal or cochain.Dual");
        }

        pub fn exteriorDerivative(
            self: *Self,
            comptime Duality: type,
            comptime k: comptime_int,
        ) *const exterior_derivative_mod.AssembledExteriorDerivative(
            cochain.Cochain(MeshType, k, Duality),
        ) {
            if (k < 0 or k >= MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no exterior derivative degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            if (Duality == cochain.Primal) {
                std.debug.assert(self.exterior_derivative_primal[k] != null);
                return &self.exterior_derivative_primal[k].?;
            }

            if (Duality == cochain.Dual) {
                std.debug.assert(self.exterior_derivative_dual[k] != null);
                return &self.exterior_derivative_dual[k].?;
            }

            @compileError("ExteriorDerivative duality must be cochain.Primal or cochain.Dual");
        }

        pub fn withHodgeStar(self: *Self, comptime k: comptime_int) !void {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no Hodge star degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            if (self.hodge_stars[k] == null) {
                self.hodge_stars[k] = try hodge_star_mod.assemble_for_degree(MeshType, k, self.allocator, self.mesh);
            }
        }

        pub fn hodgeStar(self: *Self, comptime k: comptime_int) *const hodge_star_mod.AssembledHodgeStar(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no Hodge star degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }
            std.debug.assert(self.hodge_stars[k] != null);
            return &self.hodge_stars[k].?;
        }

        pub fn withHodgeStarInverse(self: *Self, comptime k: comptime_int) !void {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no inverse Hodge star degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            if (self.hodge_star_inverses[k] == null) {
                self.hodge_star_inverses[k] = try hodge_star_mod.assemble_inverse_for_degree(MeshType, k, self.allocator, self.mesh);
            }
        }

        pub fn hodgeStarInverse(self: *Self, comptime k: comptime_int) *const hodge_star_mod.AssembledHodgeStarInverse(
            cochain.Cochain(MeshType, MeshType.topological_dimension - k, cochain.Dual),
        ) {
            if (k < 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no inverse Hodge star degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }
            std.debug.assert(self.hodge_star_inverses[k] != null);
            return &self.hodge_star_inverses[k].?;
        }

        pub fn withCodifferential(self: *Self, comptime k: comptime_int) !void {
            if (k <= 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no codifferential degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }

            const slot_index = k - 1;
            if (self.codifferentials[slot_index] == null) {
                self.codifferentials[slot_index] = try codifferential_mod.assemble_for_degree(MeshType, k, self.allocator, self.mesh);
            }
        }

        pub fn codifferential(self: *Self, comptime k: comptime_int) *const codifferential_mod.AssembledCodifferential(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            if (k <= 0 or k > MeshType.topological_dimension) {
                @compileError(std.fmt.comptimePrint(
                    "no codifferential degree {d} on a {d}-dimensional mesh",
                    .{ k, MeshType.topological_dimension },
                ));
            }
            std.debug.assert(self.codifferentials[k - 1] != null);
            return &self.codifferentials[k - 1].?;
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
