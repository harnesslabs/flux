//! Operator assembly context bound to a single mesh.
//!
//! The mesh owns topology and geometric primitives. The operator context owns
//! only the assembled operator state actually requested by a given problem.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
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
/// The public API is a heap-allocated handle so downstream structs can store a
/// pointer rather than duplicating ownership by value.
pub fn OperatorContext(comptime MeshType: type) type {
    comptime {
        @setEvalBranchQuota(20_000);
    }

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

        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !*Self {
            const self = try allocator.create(Self);
            errdefer allocator.destroy(self);

            self.* = .{
                .allocator = allocator,
                .mesh = mesh,
                .exterior_derivative_primal = std.mem.zeroes(ExteriorDerivativePrimalSlots),
                .exterior_derivative_dual = std.mem.zeroes(ExteriorDerivativeDualSlots),
                .hodge_stars = std.mem.zeroes(HodgeStarSlots),
                .hodge_star_inverses = std.mem.zeroes(HodgeStarInverseSlots),
                .codifferentials = std.mem.zeroes(CodifferentialSlots),
                .laplacians = std.mem.zeroes(LaplacianSlots),
            };
            return self;
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
                if (self.laplacians[k]) |*op| op.deinit(self.allocator);
            }
            self.allocator.destroy(self);
        }

        fn ensureExteriorDerivative(self: *Self, comptime Duality: type, comptime k: comptime_int) !void {
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
        ) !*const exterior_derivative_mod.AssembledExteriorDerivative(
            cochain.Cochain(MeshType, k, Duality),
        ) {
            try self.ensureExteriorDerivative(Duality, k);

            if (Duality == cochain.Primal) {
                return &self.exterior_derivative_primal[k].?;
            }

            if (Duality == cochain.Dual) {
                return &self.exterior_derivative_dual[k].?;
            }

            @compileError("ExteriorDerivative duality must be cochain.Primal or cochain.Dual");
        }

        fn ensureHodgeStar(self: *Self, comptime k: comptime_int) !void {
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

        pub fn hodgeStar(self: *Self, comptime k: comptime_int) !*const hodge_star_mod.AssembledHodgeStar(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            try self.ensureHodgeStar(k);
            return &self.hodge_stars[k].?;
        }

        fn ensureHodgeStarInverse(self: *Self, comptime k: comptime_int) !void {
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

        pub fn hodgeStarInverse(self: *Self, comptime k: comptime_int) !*const hodge_star_mod.AssembledHodgeStarInverse(
            cochain.Cochain(MeshType, MeshType.topological_dimension - k, cochain.Dual),
        ) {
            try self.ensureHodgeStarInverse(k);
            return &self.hodge_star_inverses[k].?;
        }

        fn ensureCodifferential(self: *Self, comptime k: comptime_int) !void {
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

        pub fn codifferential(self: *Self, comptime k: comptime_int) !*const codifferential_mod.AssembledCodifferential(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            try self.ensureCodifferential(k);
            return &self.codifferentials[k - 1].?;
        }

        fn ensureLaplacian(self: *Self, comptime k: comptime_int) !void {
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

        /// Access an assembled Δₖ, assembling it on first use.
        pub fn laplacian(self: *Self, comptime k: comptime_int) !*const laplacian_mod.AssembledLaplacian(
            cochain.Cochain(MeshType, k, cochain.Primal),
        ) {
            try self.ensureLaplacian(k);
            return &self.laplacians[k].?;
        }
    };
}

const Mesh2D = topology.Mesh(2, 2);

test "OperatorContext deinit releases assembled operator storage" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    {
        var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
        defer mesh.deinit(allocator);

        const operator_context = try OperatorContext(Mesh2D).init(allocator, &mesh);
        defer operator_context.deinit();

        _ = try operator_context.exteriorDerivative(cochain.Primal, 0);
        _ = try operator_context.exteriorDerivative(cochain.Primal, 1);
        _ = try operator_context.exteriorDerivative(cochain.Dual, 0);
        _ = try operator_context.hodgeStar(0);
        _ = try operator_context.hodgeStar(1);
        _ = try operator_context.hodgeStar(2);
        _ = try operator_context.hodgeStarInverse(1);
        _ = try operator_context.codifferential(1);
        _ = try operator_context.laplacian(0);
        _ = try operator_context.laplacian(1);
        _ = try operator_context.laplacian(2);
    }

    try testing.expect(gpa.deinit() == .ok);
}
