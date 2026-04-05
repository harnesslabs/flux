//! Composable boundary conditions for cochains on meshes with boundary.
//!
//! Each boundary condition is a value with a single public operation:
//! `apply(allocator, input) -> output`. The input and output cochain types are
//! identical, so boundary conditions compose without changing the field space.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const exterior_derivative_mod = @import("exterior_derivative.zig");

fn boundaryEffectiveDegree(comptime InputType: type) comptime_int {
    return if (InputType.duality == cochain.Dual)
        InputType.MeshT.topological_dimension - InputType.degree
    else
        InputType.degree;
}

fn BoundarySelection(comptime InputType: type) type {
    return struct {
        const Self = @This();

        mask: []bool,

        pub fn initBoundary(allocator: std.mem.Allocator, mesh: *const InputType.MeshT) !Self {
            _ = allocator;
            _ = mesh;
            @panic("not yet implemented");
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.mask);
        }

        pub fn contains(self: Self, idx: u32) bool {
            return self.mask[idx];
        }
    };
}

pub fn BoundaryCondition(comptime BC: type, comptime InputType: type) void {
    comptime {
        if (!@hasDecl(BC, "apply")) {
            @compileError("boundary condition must declare apply(allocator, input)");
        }

        const apply_type = @TypeOf(BC.apply);
        const info = @typeInfo(apply_type);
        if (info != .@"fn") {
            @compileError("boundary condition apply must be a function");
        }

        const fn_info = info.@"fn";
        if (fn_info.params.len != 3) {
            @compileError("boundary condition apply must accept self, allocator, input");
        }
        if (fn_info.return_type == null) {
            @compileError("boundary condition apply must return an error union");
        }

        const return_info = @typeInfo(fn_info.return_type.?);
        if (return_info != .error_union) {
            @compileError("boundary condition apply must return an error union");
        }
        if (return_info.error_union.payload != InputType) {
            @compileError("boundary condition apply must return the same cochain type");
        }
    }
}

pub fn Dirichlet(comptime InputType: type) type {
    return struct {
        const Self = @This();

        selection: BoundarySelection(InputType),
        boundary_values: []const f64,

        pub fn initBoundary(
            allocator: std.mem.Allocator,
            mesh: *const InputType.MeshT,
            boundary_values: []const f64,
        ) !Self {
            std.debug.assert(boundary_values.len == InputType.num_cells(mesh));
            return .{
                .selection = try BoundarySelection(InputType).initBoundary(allocator, mesh),
                .boundary_values = boundary_values,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.selection.deinit(allocator);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            _ = self;
            _ = allocator;
            _ = input;
            @panic("not yet implemented");
        }
    };
}

pub fn PEC(comptime InputType: type) type {
    return struct {
        const Self = @This();

        inner: Dirichlet(InputType),

        pub fn initBoundary(allocator: std.mem.Allocator, mesh: *const InputType.MeshT) !Self {
            const zero_values = try allocator.alloc(f64, InputType.num_cells(mesh));
            @memset(zero_values, 0.0);
            return .{
                .inner = try Dirichlet(InputType).initBoundary(allocator, mesh, zero_values),
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.inner.boundary_values);
            self.inner.deinit(allocator);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            return self.inner.apply(allocator, input);
        }
    };
}

pub fn NoSlip(comptime InputType: type) type {
    return struct {
        const Self = @This();

        inner: Dirichlet(InputType),

        pub fn initBoundary(allocator: std.mem.Allocator, mesh: *const InputType.MeshT) !Self {
            const zero_values = try allocator.alloc(f64, InputType.num_cells(mesh));
            @memset(zero_values, 0.0);
            return .{
                .inner = try Dirichlet(InputType).initBoundary(allocator, mesh, zero_values),
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.inner.boundary_values);
            self.inner.deinit(allocator);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            return self.inner.apply(allocator, input);
        }
    };
}

pub fn Periodic(comptime InputType: type) type {
    return struct {
        const Self = @This();

        representatives: []u32,

        pub fn init(allocator: std.mem.Allocator, mesh: *const InputType.MeshT, representatives: []const u32) !Self {
            std.debug.assert(representatives.len == InputType.num_cells(mesh));
            const owned = try allocator.dupe(u32, representatives);
            return .{ .representatives = owned };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.representatives);
        }

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: InputType) !InputType {
            _ = self;
            _ = allocator;
            _ = input;
            @panic("not yet implemented");
        }
    };
}

pub fn compose(first: anytype, second: anytype) type {
    const First = @TypeOf(first);
    const Second = @TypeOf(second);

    return struct {
        const Self = @This();

        first: First,
        second: Second,

        pub fn apply(self: Self, allocator: std.mem.Allocator, input: anytype) !@TypeOf(input) {
            var intermediate = try self.first.apply(allocator, input);
            defer intermediate.deinit(allocator);
            return try self.second.apply(allocator, intermediate);
        }
    };
}

const Mesh2D = topology.Mesh(2, 2);
const Mesh3D = topology.Mesh(3, 3);
const Primal0_2D = cochain.Cochain(Mesh2D, 0, cochain.Primal);
const Primal1_2D = cochain.Cochain(Mesh2D, 1, cochain.Primal);
const Primal1_3D = cochain.Cochain(Mesh3D, 1, cochain.Primal);

test "boundary condition concept accepts Dirichlet and periodic on cochains" {
    BoundaryCondition(Dirichlet(Primal0_2D), Primal0_2D);
    BoundaryCondition(PEC(Primal1_2D), Primal1_2D);
    BoundaryCondition(NoSlip(Primal1_3D), Primal1_3D);
    BoundaryCondition(Periodic(Primal0_2D), Primal0_2D);
}

test "boundary selection marks boundary vertices on 2D mesh" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var selection = try BoundarySelection(Primal0_2D).initBoundary(allocator, &mesh);
    defer selection.deinit(allocator);

    try testing.expect(selection.contains(0));
    try testing.expect(selection.contains(1));
    try testing.expect(selection.contains(2));
    try testing.expect(selection.contains(3));
    try testing.expect(!selection.contains(4));
    try testing.expect(selection.contains(5));
    try testing.expect(selection.contains(6));
    try testing.expect(selection.contains(7));
    try testing.expect(selection.contains(8));
}

test "boundary selection marks boundary edges on 3D mesh" {
    const allocator = testing.allocator;
    var mesh = try Mesh3D.uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var selection = try BoundarySelection(Primal1_3D).initBoundary(allocator, &mesh);
    defer selection.deinit(allocator);

    var selected_count: usize = 0;
    for (selection.mask) |selected| {
        if (selected) selected_count += 1;
    }

    try testing.expect(selected_count > 0);
    try testing.expect(selected_count < mesh.num_edges());
}

test "Dirichlet is idempotent on boundary vertices" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var field = try Primal0_2D.init(allocator, &mesh);
    defer field.deinit(allocator);
    for (field.values, 0..) |*value, idx| {
        value.* = @floatFromInt(idx);
    }

    const boundary_values = try allocator.alloc(f64, field.values.len);
    defer allocator.free(boundary_values);
    for (boundary_values, 0..) |*value, idx| {
        value.* = 100.0 + @as(f64, @floatFromInt(idx));
    }

    var bc = try Dirichlet(Primal0_2D).initBoundary(allocator, &mesh, boundary_values);
    defer bc.deinit(allocator);

    var once = try bc.apply(allocator, field);
    defer once.deinit(allocator);

    var twice = try bc.apply(allocator, once);
    defer twice.deinit(allocator);

    for (once.values, twice.values) |lhs, rhs| {
        try testing.expectEqual(lhs, rhs);
    }
}

test "PEC is idempotent on boundary edges" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var field = try Primal1_2D.init(allocator, &mesh);
    defer field.deinit(allocator);
    for (field.values, 0..) |*value, idx| {
        value.* = 1.0 + @as(f64, @floatFromInt(idx));
    }

    var bc = try PEC(Primal1_2D).initBoundary(allocator, &mesh);
    defer bc.deinit(allocator);

    var once = try bc.apply(allocator, field);
    defer once.deinit(allocator);

    var twice = try bc.apply(allocator, once);
    defer twice.deinit(allocator);

    for (once.values, twice.values) |lhs, rhs| {
        try testing.expectEqual(lhs, rhs);
    }
}

test "periodic averages paired boundary vertices and is idempotent" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var field = try Primal0_2D.init(allocator, &mesh);
    defer field.deinit(allocator);
    field.values[0] = 2.0;
    field.values[1] = 4.0;
    field.values[2] = 10.0;
    field.values[3] = 14.0;

    const representatives = [_]u32{ 0, 1, 0, 1 };
    var bc = try Periodic(Primal0_2D).init(allocator, &mesh, &representatives);
    defer bc.deinit(allocator);

    var once = try bc.apply(allocator, field);
    defer once.deinit(allocator);

    var twice = try bc.apply(allocator, once);
    defer twice.deinit(allocator);

    try testing.expectEqual(@as(f64, 6.0), once.values[0]);
    try testing.expectEqual(@as(f64, 6.0), once.values[2]);
    try testing.expectEqual(@as(f64, 9.0), once.values[1]);
    try testing.expectEqual(@as(f64, 9.0), once.values[3]);

    for (once.values, twice.values) |lhs, rhs| {
        try testing.expectEqual(lhs, rhs);
    }
}

test "Dirichlet composed with exterior derivative zeros boundary-edge gradient for constant boundary data" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var field = try Primal0_2D.init(allocator, &mesh);
    defer field.deinit(allocator);

    const coords = mesh.vertices.slice().items(.coords);
    for (field.values, coords) |*value, coord| {
        value.* = coord[0] + 2.0 * coord[1];
    }

    const boundary_values = try allocator.alloc(f64, field.values.len);
    defer allocator.free(boundary_values);
    @memset(boundary_values, 5.0);

    var bc = try Dirichlet(Primal0_2D).initBoundary(allocator, &mesh, boundary_values);
    defer bc.deinit(allocator);

    var constrained = try bc.apply(allocator, field);
    defer constrained.deinit(allocator);

    var derivative = try exterior_derivative_mod.exterior_derivative(allocator, constrained);
    defer derivative.deinit(allocator);

    for (mesh.boundary_edges) |edge_idx| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), derivative.values[edge_idx], 1e-12);
    }
}
