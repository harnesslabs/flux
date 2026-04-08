//! VTK Unstructured Grid (`.vtu`) serializer.
//!
//! Writes mesh geometry and attached field data (0-forms as PointData, higher
//! forms as CellData) to the VTK XML format for inspection in ParaView. All
//! numeric data is ASCII — portable and human-readable, not optimized for size.
//!
//! Reference: VTK File Formats specification, §19.3 "XML File Formats"
//! (Kitware, The VTK User's Guide, 11th Edition).

const std = @import("std");
const testing = std.testing;
const topology = @import("../topology/mesh.zig");

/// VTK cell type for triangles (linear triangle, type 5 in the VTK spec).
const vtk_triangle: u8 = 5;
/// VTK cell type for tetrahedra (linear tetrahedron, type 10 in the VTK spec).
const vtk_tetrahedron: u8 = 10;

/// Named field to attach to VTK output as PointData or CellData.
pub const DataArraySlice = struct {
    name: []const u8,
    values: []const f64,
    num_components: u8 = 1,
};

/// Write a VTK Unstructured Grid (.vtu) XML file for a simplicial mesh.
///
/// The .vtu format stores:
/// - `<Points>`: vertex coordinates (3-component, padded with zeros for embeddings < 3)
/// - `<Cells>`: top-dimensional simplex connectivity, offsets, and VTK types
/// - `<PointData>`: scalar or vector fields on vertices
/// - `<CellData>`: scalar or vector fields on top-dimensional cells
///
/// All numeric data is written in ASCII. This is not the fastest format
/// but is portable, human-readable, and sufficient for visualization.
///
/// Reference: VTK File Formats specification, §19.3 "XML File Formats"
/// (Kitware, The VTK User's Guide, 11th Edition).
pub fn write(
    writer: anytype,
    mesh: anytype,
    point_data: []const DataArraySlice,
    cell_data: []const DataArraySlice,
) !void {
    const MeshType = @TypeOf(mesh);
    comptime {
        if (!@hasDecl(MeshType, "embedding_dimension") or !@hasDecl(MeshType, "topological_dimension")) {
            @compileError("vtk.write requires a mesh type with embedding_dimension and topological_dimension");
        }
    }

    const embedding_dimension = MeshType.embedding_dimension;
    const topological_dimension = MeshType.topological_dimension;
    const num_vertices = mesh.num_vertices();
    const num_cells = cellCount(mesh, topological_dimension);

    // -- XML header and VTKFile element --
    try writer.writeAll("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    try writer.writeAll("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    try writer.writeAll("  <UnstructuredGrid>\n");
    try writer.print("    <Piece NumberOfPoints=\"{d}\" NumberOfCells=\"{d}\">\n", .{ num_vertices, num_cells });

    // -- PointData (0-forms: one value per vertex) --
    if (point_data.len > 0) {
        try writer.writeAll("      <PointData>\n");
        for (point_data) |pd| {
            std.debug.assert(pd.num_components > 0);
            std.debug.assert(pd.values.len == num_vertices * pd.num_components);
            try writeDataArray(writer, pd.name, pd.values, pd.num_components);
        }
        try writer.writeAll("      </PointData>\n");
    }

    // -- CellData (1-forms or 2-forms mapped to faces) --
    if (cell_data.len > 0) {
        try writer.writeAll("      <CellData>\n");
        for (cell_data) |cd| {
            std.debug.assert(cd.num_components > 0);
            std.debug.assert(cd.values.len == num_cells * cd.num_components);
            try writeDataArray(writer, cd.name, cd.values, cd.num_components);
        }
        try writer.writeAll("      </CellData>\n");
    }

    // -- Points --
    try writer.writeAll("      <Points>\n");
    try writer.writeAll("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    {
        const coords = mesh.vertices.slice().items(.coords);
        for (coords) |coord| {
            try writer.writeAll("          ");
            // Always write 3 components; pad with zeros for dimensions < 3.
            inline for (0..3) |d| {
                if (d > 0) try writer.writeByte(' ');
                if (d < embedding_dimension) {
                    try writer.print("{e}", .{coord[d]});
                } else {
                    try writer.writeAll("0e0");
                }
            }
            try writer.writeByte('\n');
        }
    }
    try writer.writeAll("        </DataArray>\n");
    try writer.writeAll("      </Points>\n");

    // -- Cells --
    try writer.writeAll("      <Cells>\n");

    // Connectivity: vertex indices for each top-dimensional cell.
    try writer.writeAll("        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");
    try writeConnectivity(writer, mesh, topological_dimension);
    try writer.writeAll("        </DataArray>\n");

    // Offsets: cumulative vertex count per cell.
    try writer.writeAll("        <DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
    try writeOffsets(writer, num_cells, topological_dimension);
    try writer.writeAll("        </DataArray>\n");

    // Types: one VTK cell type code per top-dimensional cell.
    try writer.writeAll("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    try writeCellTypes(writer, num_cells, topological_dimension);
    try writer.writeAll("        </DataArray>\n");

    try writer.writeAll("      </Cells>\n");

    // -- Close elements --
    try writer.writeAll("    </Piece>\n");
    try writer.writeAll("  </UnstructuredGrid>\n");
    try writer.writeAll("</VTKFile>\n");
}

fn cellCount(mesh: anytype, comptime topological_dimension: usize) u32 {
    return switch (topological_dimension) {
        2 => mesh.num_faces(),
        3 => mesh.num_tets(),
        else => @compileError("VTK export currently supports only topological dimensions 2 and 3"),
    };
}

fn verticesPerCell(comptime topological_dimension: usize) u32 {
    return switch (topological_dimension) {
        2 => 3,
        3 => 4,
        else => @compileError("VTK export currently supports only topological dimensions 2 and 3"),
    };
}

fn cellTypeForDimension(comptime topological_dimension: usize) u8 {
    return switch (topological_dimension) {
        2 => vtk_triangle,
        3 => vtk_tetrahedron,
        else => @compileError("VTK export currently supports only topological dimensions 2 and 3"),
    };
}

fn writeConnectivity(writer: anytype, mesh: anytype, comptime topological_dimension: usize) !void {
    switch (topological_dimension) {
        2 => {
            const face_verts = mesh.simplices(2).items(.vertices);
            for (face_verts) |verts| {
                try writer.print("          {d} {d} {d}\n", .{ verts[0], verts[1], verts[2] });
            }
        },
        3 => {
            const tet_verts = mesh.simplices(3).items(.vertices);
            for (tet_verts) |verts| {
                try writer.print("          {d} {d} {d} {d}\n", .{ verts[0], verts[1], verts[2], verts[3] });
            }
        },
        else => @compileError("VTK export currently supports only topological dimensions 2 and 3"),
    }
}

fn writeOffsets(writer: anytype, num_cells: u32, comptime topological_dimension: usize) !void {
    try writer.writeAll("          ");
    for (0..num_cells) |cell_index| {
        if (cell_index > 0) try writer.writeByte(' ');
        try writer.print("{d}", .{(cell_index + 1) * verticesPerCell(topological_dimension)});
    }
    try writer.writeByte('\n');
}

fn writeCellTypes(writer: anytype, num_cells: u32, comptime topological_dimension: usize) !void {
    try writer.writeAll("          ");
    for (0..num_cells) |cell_index| {
        if (cell_index > 0) try writer.writeByte(' ');
        try writer.print("{d}", .{cellTypeForDimension(topological_dimension)});
    }
    try writer.writeByte('\n');
}

/// Write a single `<DataArray>` element with Float64 scalar or vector data.
fn writeDataArray(writer: anytype, name: []const u8, values: []const f64, num_components: u8) !void {
    std.debug.assert(num_components > 0);
    if (num_components == 1) {
        try writer.print("        <DataArray type=\"Float64\" Name=\"{s}\" format=\"ascii\">\n", .{name});
    } else {
        try writer.print(
            "        <DataArray type=\"Float64\" Name=\"{s}\" NumberOfComponents=\"{d}\" format=\"ascii\">\n",
            .{ name, num_components },
        );
    }
    try writer.writeAll("          ");
    for (values, 0..) |v, i| {
        if (i > 0) try writer.writeByte(' ');
        try writer.print("{e}", .{v});
    }
    try writer.writeByte('\n');
    try writer.writeAll("        </DataArray>\n");
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

test "vtu output contains valid XML structure for minimal mesh" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const xml = output.items;

    // Must contain XML declaration
    try testing.expect(std.mem.startsWith(u8, xml, "<?xml version=\"1.0\""));
    // Must contain VTKFile root element
    try testing.expect(std.mem.indexOf(u8, xml, "<VTKFile type=\"UnstructuredGrid\"") != null);
    // Must contain correct point/cell counts: 4 vertices, 2 faces
    try testing.expect(std.mem.indexOf(u8, xml, "NumberOfPoints=\"4\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "NumberOfCells=\"2\"") != null);
    // Must contain Points, Cells, connectivity, offsets, types
    try testing.expect(std.mem.indexOf(u8, xml, "<Points>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "<Cells>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"connectivity\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"offsets\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"types\"") != null);
    // Must close all elements
    try testing.expect(std.mem.indexOf(u8, xml, "</VTKFile>") != null);
}

test "vtu vertex coordinates round-trip to machine precision" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 2, 2, 3.0, 5.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    // Parse back all vertex coordinates from the Points DataArray.
    // The coordinates appear between <DataArray ...Float64...NumberOfComponents="3"...>
    // and the closing </DataArray> inside <Points>.
    const xml = output.items;
    const points_start_tag = "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    const points_start = std.mem.indexOf(u8, xml, points_start_tag).? + points_start_tag.len;
    const points_end = std.mem.indexOfPos(u8, xml, points_start, "</DataArray>").?;
    const points_block = xml[points_start..points_end];

    const coords = mesh.vertices.slice().items(.coords);
    var vertex_idx: usize = 0;
    var line_iter = std.mem.tokenizeScalar(u8, points_block, '\n');
    while (line_iter.next()) |line| {
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        var component_iter = std.mem.tokenizeScalar(u8, trimmed, ' ');
        // x
        const x_str = component_iter.next().?;
        const x = try std.fmt.parseFloat(f64, x_str);
        try testing.expectApproxEqAbs(coords[vertex_idx][0], x, 1e-15);
        // y
        const y_str = component_iter.next().?;
        const y = try std.fmt.parseFloat(f64, y_str);
        try testing.expectApproxEqAbs(coords[vertex_idx][1], y, 1e-15);
        // z (should be 0 for 2D)
        const z_str = component_iter.next().?;
        const z = try std.fmt.parseFloat(f64, z_str);
        try testing.expectApproxEqAbs(@as(f64, 0.0), z, 1e-15);

        vertex_idx += 1;
    }
    try testing.expectEqual(mesh.num_vertices(), @as(u32, @intCast(vertex_idx)));
}

test "vtu 0-form PointData round-trips to machine precision" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // Create a 0-form: f(v) = x² + y at each vertex.
    const coords = mesh.vertices.slice().items(.coords);
    const values = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(values);
    for (coords, 0..) |c, i| {
        values[i] = c[0] * c[0] + c[1];
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const pd = [_]DataArraySlice{.{ .name = "temperature", .values = values }};
    try write(output.writer(allocator), mesh, &pd, &.{});

    const xml = output.items;

    // Verify PointData section exists
    try testing.expect(std.mem.indexOf(u8, xml, "<PointData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"temperature\"") != null);

    // Parse values back and compare
    const parsed = try parseDataArray(allocator, xml, "temperature");
    defer allocator.free(parsed);

    try testing.expectEqual(values.len, parsed.len);
    for (values, parsed) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

test "vtu 2-form CellData round-trips to machine precision" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // Create a 2-form: assign face index as value (distinct values for verification).
    const values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(values);
    for (values, 0..) |*v, i| {
        v.* = @as(f64, @floatFromInt(i)) * 1.5 + 0.7;
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const cd = [_]DataArraySlice{.{ .name = "vorticity", .values = values }};
    try write(output.writer(allocator), mesh, &.{}, &cd);

    const xml = output.items;

    // Verify CellData section exists
    try testing.expect(std.mem.indexOf(u8, xml, "<CellData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"vorticity\"") != null);

    // Parse values back and compare
    const parsed = try parseDataArray(allocator, xml, "vorticity");
    defer allocator.free(parsed);

    try testing.expectEqual(values.len, parsed.len);
    for (values, parsed) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

test "vtu with both PointData and CellData" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const point_values = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(point_values);
    for (point_values, 0..) |*v, i| {
        v.* = @as(f64, @floatFromInt(i));
    }

    const cell_values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(cell_values);
    for (cell_values, 0..) |*v, i| {
        v.* = @as(f64, @floatFromInt(i)) * 10.0;
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const pd = [_]DataArraySlice{.{ .name = "pressure", .values = point_values }};
    const cd = [_]DataArraySlice{.{ .name = "flux", .values = cell_values }};
    try write(output.writer(allocator), mesh, &pd, &cd);

    const xml = output.items;
    try testing.expect(std.mem.indexOf(u8, xml, "<PointData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "<CellData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"pressure\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"flux\"") != null);
}

test "vtu triangle connectivity matches mesh face vertices" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const xml = output.items;

    // Parse connectivity array
    const conn_start_tag = "Name=\"connectivity\" format=\"ascii\">";
    const conn_start = std.mem.indexOf(u8, xml, conn_start_tag).? + conn_start_tag.len;
    const conn_end = std.mem.indexOfPos(u8, xml, conn_start, "</DataArray>").?;
    const conn_block = xml[conn_start..conn_end];

    const face_verts = mesh.simplices(2).items(.vertices);
    var face_idx: usize = 0;
    var line_iter = std.mem.tokenizeScalar(u8, conn_block, '\n');
    while (line_iter.next()) |line| {
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        var tok = std.mem.tokenizeScalar(u8, trimmed, ' ');
        const v0 = try std.fmt.parseInt(u32, tok.next().?, 10);
        const v1 = try std.fmt.parseInt(u32, tok.next().?, 10);
        const v2 = try std.fmt.parseInt(u32, tok.next().?, 10);

        try testing.expectEqual(face_verts[face_idx][0], v0);
        try testing.expectEqual(face_verts[face_idx][1], v1);
        try testing.expectEqual(face_verts[face_idx][2], v2);
        face_idx += 1;
    }
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(face_idx)));
}

test "vtu all cell types are triangle (type 5)" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 3, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const xml = output.items;

    // Parse types array
    const types_start_tag = "Name=\"types\" format=\"ascii\">";
    const types_start = std.mem.indexOf(u8, xml, types_start_tag).? + types_start_tag.len;
    const types_end = std.mem.indexOfPos(u8, xml, types_start, "</DataArray>").?;
    const types_block = xml[types_start..types_end];

    var count: usize = 0;
    var tok = std.mem.tokenizeAny(u8, types_block, " \t\n\r");
    while (tok.next()) |t| {
        const val = try std.fmt.parseInt(u8, t, 10);
        try testing.expectEqual(vtk_triangle, val);
        count += 1;
    }
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(count)));
}

test "vtu 3D piece counts use tetrahedra as cells" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const expected_piece = try std.fmt.allocPrint(allocator, "NumberOfPoints=\"{d}\" NumberOfCells=\"{d}\"", .{
        mesh.num_vertices(),
        mesh.num_tets(),
    });
    defer allocator.free(expected_piece);

    try testing.expect(std.mem.indexOf(u8, output.items, expected_piece) != null);
}

test "vtu tetrahedral connectivity matches mesh tet vertices" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const xml = output.items;
    const conn_start_tag = "Name=\"connectivity\" format=\"ascii\">";
    const conn_start = std.mem.indexOf(u8, xml, conn_start_tag).? + conn_start_tag.len;
    const conn_end = std.mem.indexOfPos(u8, xml, conn_start, "</DataArray>").?;
    const conn_block = xml[conn_start..conn_end];

    const tet_verts = mesh.simplices(3).items(.vertices);
    var tet_idx: usize = 0;
    var line_iter = std.mem.tokenizeScalar(u8, conn_block, '\n');
    while (line_iter.next()) |line| {
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        var tok = std.mem.tokenizeScalar(u8, trimmed, ' ');
        const token0 = tok.next();
        const token1 = tok.next();
        const token2 = tok.next();
        const token3 = tok.next();
        try testing.expect(token0 != null);
        try testing.expect(token1 != null);
        try testing.expect(token2 != null);
        try testing.expect(token3 != null);

        const v0 = try std.fmt.parseInt(u32, token0.?, 10);
        const v1 = try std.fmt.parseInt(u32, token1.?, 10);
        const v2 = try std.fmt.parseInt(u32, token2.?, 10);
        const v3 = try std.fmt.parseInt(u32, token3.?, 10);

        try testing.expectEqual(tet_verts[tet_idx][0], v0);
        try testing.expectEqual(tet_verts[tet_idx][1], v1);
        try testing.expectEqual(tet_verts[tet_idx][2], v2);
        try testing.expectEqual(tet_verts[tet_idx][3], v3);
        try testing.expect(tok.next() == null);
        tet_idx += 1;
    }

    try testing.expectEqual(mesh.num_tets(), @as(u32, @intCast(tet_idx)));
}

test "vtu 3D cell types are tetrahedron (type 10)" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write(output.writer(allocator), mesh, &.{}, &.{});

    const xml = output.items;
    const types_start_tag = "Name=\"types\" format=\"ascii\">";
    const types_start = std.mem.indexOf(u8, xml, types_start_tag).? + types_start_tag.len;
    const types_end = std.mem.indexOfPos(u8, xml, types_start, "</DataArray>").?;
    const types_block = xml[types_start..types_end];

    var count: usize = 0;
    var tok = std.mem.tokenizeAny(u8, types_block, " \t\n\r");
    while (tok.next()) |t| {
        const val = try std.fmt.parseInt(u8, t, 10);
        try testing.expectEqual(@as(u8, 10), val);
        count += 1;
    }

    try testing.expectEqual(mesh.num_tets(), @as(u32, @intCast(count)));
}

test "vtu 3D PointData and CellData round-trip with tetrahedral mesh" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const point_values = try allocator.alloc(f64, mesh.num_vertices());
    defer allocator.free(point_values);
    for (point_values, 0..) |*value, i| {
        value.* = @as(f64, @floatFromInt(i)) + 0.25;
    }

    const cell_values = try allocator.alloc(f64, mesh.num_tets());
    defer allocator.free(cell_values);
    for (cell_values, 0..) |*value, i| {
        value.* = @as(f64, @floatFromInt(i)) * 2.0 + 0.5;
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const pd = [_]DataArraySlice{.{ .name = "temperature_3d", .values = point_values }};
    const cd = [_]DataArraySlice{.{ .name = "pressure_3d", .values = cell_values }};
    try write(output.writer(allocator), mesh, &pd, &cd);

    const parsed_point = try parseDataArray(allocator, output.items, "temperature_3d");
    defer allocator.free(parsed_point);
    try testing.expectEqual(point_values.len, parsed_point.len);
    for (point_values, parsed_point) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }

    const parsed_cell = try parseDataArray(allocator, output.items, "pressure_3d");
    defer allocator.free(parsed_cell);
    try testing.expectEqual(cell_values.len, parsed_cell.len);
    for (cell_values, parsed_cell) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

test "vtu 3D vector PointData and CellData preserve component count and layout" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(3, 3).uniform_tetrahedral_grid(allocator, 1, 1, 1, 1.0, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const point_values = try allocator.alloc(f64, mesh.num_vertices() * 3);
    defer allocator.free(point_values);
    for (point_values, 0..) |*value, i| {
        value.* = @as(f64, @floatFromInt(i)) * 0.25;
    }

    const cell_values = try allocator.alloc(f64, mesh.num_tets() * 3);
    defer allocator.free(cell_values);
    for (cell_values, 0..) |*value, i| {
        value.* = @as(f64, @floatFromInt(i)) * 0.5;
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const pd = [_]DataArraySlice{.{ .name = "velocity_3d", .values = point_values, .num_components = 3 }};
    const cd = [_]DataArraySlice{.{ .name = "flux_3d", .values = cell_values, .num_components = 3 }};
    try write(output.writer(allocator), mesh, &pd, &cd);

    try testing.expect(std.mem.indexOf(u8, output.items, "Name=\"velocity_3d\" NumberOfComponents=\"3\"") != null);
    try testing.expect(std.mem.indexOf(u8, output.items, "Name=\"flux_3d\" NumberOfComponents=\"3\"") != null);

    const parsed_point = try parseDataArray(allocator, output.items, "velocity_3d");
    defer allocator.free(parsed_point);
    try testing.expectEqual(point_values.len, parsed_point.len);
    for (point_values, parsed_point) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }

    const parsed_cell = try parseDataArray(allocator, output.items, "flux_3d");
    defer allocator.free(parsed_cell);
    try testing.expectEqual(cell_values.len, parsed_cell.len);
    for (cell_values, parsed_cell) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Parsing
// ───────────────────────────────────────────────────────────────────────────

/// Check whether an XML tag fragment contains `Name="<expected>"`.
/// Avoids heap allocation by scanning for the prefix and comparing
/// the name bytes inline.
fn matchNameAttribute(tag: []const u8, expected: []const u8) bool {
    const prefix = "Name=\"";
    var pos: usize = 0;
    while (std.mem.indexOfPos(u8, tag, pos, prefix)) |attr_start| {
        const value_start = attr_start + prefix.len;
        const value_end = value_start + expected.len;
        if (value_end < tag.len and
            std.mem.eql(u8, tag[value_start..value_end], expected) and
            tag[value_end] == '"')
        {
            return true;
        }
        pos = attr_start + 1;
    }
    return false;
}

/// Parse a named `<DataArray>` from VTK XML and return its f64 values.
///
/// Searches the XML for a `<DataArray ... Name="name" ...>` element and
/// parses its whitespace-separated body as Float64 values. Useful for
/// reading back .vtu files for verification or as a building block for
/// a full VTK reader.
///
/// Caller owns the returned slice.
pub fn parseDataArray(allocator: std.mem.Allocator, xml: []const u8, name: []const u8) ![]f64 {
    // Find the DataArray with this name
    var search_pos: usize = 0;
    while (std.mem.indexOfPos(u8, xml, search_pos, "<DataArray")) |tag_start| {
        const tag_end = std.mem.indexOfPos(u8, xml, tag_start, ">") orelse break;
        const tag = xml[tag_start..tag_end];

        // Check if this tag contains our name — match without allocating.
        // Look for Name="<name>" by finding the prefix, then comparing
        // the name bytes directly, then checking for the closing quote.
        if (matchNameAttribute(tag, name)) {
            const data_start = tag_end + 1;
            const data_end = std.mem.indexOfPos(u8, xml, data_start, "</DataArray>") orelse break;
            const data_block = xml[data_start..data_end];

            // Count values first
            var count: usize = 0;
            var counter = std.mem.tokenizeAny(u8, data_block, " \t\n\r");
            while (counter.next()) |_| count += 1;

            // Parse values
            const result = try allocator.alloc(f64, count);
            var idx: usize = 0;
            var parser = std.mem.tokenizeAny(u8, data_block, " \t\n\r");
            while (parser.next()) |token| {
                result[idx] = try std.fmt.parseFloat(f64, token);
                idx += 1;
            }
            return result;
        }

        search_pos = tag_end;
    }
    return error.DataArrayNotFound;
}

test "vtu 0-form PointData round-trips with random data across grid sizes" {
    // Property test: for multiple grid sizes, random 0-form values
    // survive write→parse round-trip to machine precision.
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0xF7C_27D_00);

    const sizes = [_][2]u32{
        .{ 1, 1 }, .{ 2, 3 }, .{ 4, 4 }, .{ 7, 5 }, .{ 10, 8 },
    };

    for (sizes) |size| {
        var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, size[0], size[1], 3.0, 2.0);
        defer mesh.deinit(allocator);

        const values = try allocator.alloc(f64, mesh.num_vertices());
        defer allocator.free(values);
        for (values) |*v| v.* = rng.random().float(f64) * 200.0 - 100.0;

        var output = std.ArrayListUnmanaged(u8){};
        defer output.deinit(allocator);

        const pd = [_]DataArraySlice{.{ .name = "random_field", .values = values }};
        try write(output.writer(allocator), mesh, &pd, &.{});

        const parsed = try parseDataArray(allocator, output.items, "random_field");
        defer allocator.free(parsed);

        try testing.expectEqual(values.len, parsed.len);
        for (values, parsed) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-15);
        }
    }
}

test "vtu 2-form CellData round-trips with random data across grid sizes" {
    // Property test: for multiple grid sizes, random 2-form values
    // survive write→parse round-trip to machine precision.
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0xF7C_27D_01);

    const sizes = [_][2]u32{
        .{ 1, 1 }, .{ 3, 2 }, .{ 5, 5 }, .{ 8, 6 },
    };

    for (sizes) |size| {
        var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, size[0], size[1], 2.0, 1.0);
        defer mesh.deinit(allocator);

        const values = try allocator.alloc(f64, mesh.num_faces());
        defer allocator.free(values);
        for (values) |*v| v.* = rng.random().float(f64) * 1000.0 - 500.0;

        var output = std.ArrayListUnmanaged(u8){};
        defer output.deinit(allocator);

        const cd = [_]DataArraySlice{.{ .name = "random_flux", .values = values }};
        try write(output.writer(allocator), mesh, &.{}, &cd);

        const parsed = try parseDataArray(allocator, output.items, "random_flux");
        defer allocator.free(parsed);

        try testing.expectEqual(values.len, parsed.len);
        for (values, parsed) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-15);
        }
    }
}
