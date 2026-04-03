//! Discrete codifferential operator δ = ★⁻¹ ∘ d ∘ ★ on primal cochains.

const std = @import("std");
const cochain = @import("../forms/cochain.zig");
const exterior_derivative = @import("exterior_derivative.zig");
const hodge_star = @import("hodge_star.zig");

fn validatePrimalInput(comptime InputType: type) void {
    if (!@hasDecl(InputType, "duality")) {
        @compileError("AssembledCodifferential requires a Cochain type");
    }
    if (InputType.duality != cochain.Primal) {
        @compileError("codifferential expects a primal cochain");
    }
    if (InputType.degree == 0) {
        @compileError("codifferential is undefined on primal 0-forms");
    }
}

fn CodifferentialResult(comptime InputType: type) type {
    return cochain.Cochain(InputType.MeshT, InputType.degree - 1, cochain.Primal);
}

pub fn AssembledCodifferential(comptime InputType: type) type {
    comptime validatePrimalInput(InputType);

    return struct {
        const Self = @This();

        mesh: *const InputType.MeshT,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            _ = self;
            _ = allocator;
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !CodifferentialResult(InputType) {
            std.debug.assert(input.mesh == self.mesh);

            var starred = try hodge_star.apply_hodge_star(allocator, input);
            defer starred.deinit(allocator);

            var dual_derivative = try exterior_derivative.apply_exterior_derivative(allocator, starred);
            defer dual_derivative.deinit(allocator);

            return hodge_star.apply_hodge_star_inverse(allocator, dual_derivative);
        }
    };
}

pub fn assemble_for_degree(
    comptime MeshType: type,
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
) !AssembledCodifferential(cochain.Cochain(MeshType, k, cochain.Primal)) {
    _ = allocator;
    return .{ .mesh = mesh };
}
