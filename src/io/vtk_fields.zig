const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");
const vtk = @import("vtk.zig");

/// Project a 1-form (one value per edge) to a per-face scalar by averaging
/// the absolute values of the three edge coefficients on each face.
pub fn project_edges_to_faces(
    allocator: std.mem.Allocator,
    mesh: anytype,
    edge_values: []const f64,
) ![]f64 {
    const MeshType = @TypeOf(mesh);
    comptime {
        if (!@hasDecl(MeshType, "topological_dimension")) {
            @compileError("project_edges_to_faces requires a mesh type with topological_dimension");
        }
        if (MeshType.topological_dimension != 2) {
            @compileError("project_edges_to_faces currently supports only 2D meshes");
        }
    }

    const num_faces = mesh.num_faces();
    const num_edges = mesh.num_edges();
    std.debug.assert(edge_values.len == num_edges);

    const face_values = try allocator.alloc(f64, num_faces);
    const boundary_2 = mesh.boundary(2);
    for (0..num_faces) |f| {
        const face_boundary = boundary_2.row(@intCast(f));
        const edge_indices = face_boundary.cols;

        var sum: f64 = 0;
        for (edge_indices) |edge_idx| {
            sum += @abs(edge_values[edge_idx]);
        }
        face_values[f] = sum / @as(f64, @floatFromInt(edge_indices.len));
    }

    return face_values;
}

/// Write a VTK snapshot of electromagnetic fields (E and B) from a Maxwell state.
pub fn write_fields(
    allocator: std.mem.Allocator,
    writer: anytype,
    mesh: anytype,
    e_values: []const f64,
    b_values: []const f64,
) !void {
    const e_projected = try project_edges_to_faces(allocator, mesh, e_values);
    defer allocator.free(e_projected);

    const cell_data = [_]vtk.DataArraySlice{
        .{ .name = "E_intensity", .values = e_projected },
        .{ .name = "B_flux", .values = b_values },
    };
    try vtk.write(writer, mesh, &.{}, &cell_data);
}

test "project_edges_to_faces averages absolute edge values per face" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).plane(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const edge_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(edge_values);
    for (edge_values, 0..) |*v, i| {
        v.* = @as(f64, @floatFromInt(i)) + 1.0;
    }

    const face_values = try project_edges_to_faces(allocator, mesh, edge_values);
    defer allocator.free(face_values);

    const boundary_2 = mesh.boundary(2);
    for (0..mesh.num_faces()) |f| {
        const face_boundary = boundary_2.row(@intCast(f));
        const edge_indices = face_boundary.cols;
        try testing.expectEqual(@as(usize, 3), edge_indices.len);

        var expected: f64 = 0;
        for (edge_indices) |e| {
            expected += @abs(edge_values[e]);
        }
        expected /= 3.0;
        try testing.expectApproxEqAbs(expected, face_values[f], 1e-15);
    }
}

test "project_edges_to_faces takes absolute values" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).plane(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const edge_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(edge_values);
    for (edge_values) |*v| v.* = -3.0;

    const face_values = try project_edges_to_faces(allocator, mesh, edge_values);
    defer allocator.free(face_values);

    for (face_values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 3.0), v, 1e-15);
    }
}

test "write_fields produces vtu with E_intensity and B_flux CellData" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).plane(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const e_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(e_values);
    for (e_values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) * 0.1;

    const b_values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(b_values);
    for (b_values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) * 0.5;

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write_fields(allocator, output.writer(allocator), mesh, e_values, b_values);

    const xml = output.items;
    try testing.expect(std.mem.indexOf(u8, xml, "<CellData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"E_intensity\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"B_flux\"") != null);

    const parsed_b = try vtk.parseDataArray(allocator, xml, "B_flux");
    defer allocator.free(parsed_b);
    try testing.expectEqual(b_values.len, parsed_b.len);
    for (b_values, parsed_b) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }

    const parsed_e = try vtk.parseDataArray(allocator, xml, "E_intensity");
    defer allocator.free(parsed_e);
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(parsed_e.len)));
}
