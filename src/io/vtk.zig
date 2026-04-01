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

/// Named field to attach to VTK output as PointData or CellData.
pub const DataArraySlice = struct {
    name: []const u8,
    values: []const f64,
};

/// Write a VTK Unstructured Grid (.vtu) XML file for a 2D mesh.
///
/// The .vtu format stores:
/// - `<Points>`: vertex coordinates (3-component, z = 0 for 2D)
/// - `<Cells>`: triangle connectivity, offsets, types
/// - `<PointData>`: scalar fields on vertices (0-forms)
/// - `<CellData>`: scalar fields on cells (1-forms on edges, 2-forms on faces)
///
/// All numeric data is written in ASCII. This is not the fastest format
/// but is portable, human-readable, and sufficient for visualization.
///
/// Reference: VTK File Formats specification, §19.3 "XML File Formats"
/// (Kitware, The VTK User's Guide, 11th Edition).
pub fn write(
    writer: anytype,
    comptime embedding_dimension: usize,
    comptime topological_dimension: usize,
    mesh: topology.Mesh(embedding_dimension, topological_dimension),
    point_data: []const DataArraySlice,
    cell_data: []const DataArraySlice,
) !void {
    const num_vertices = mesh.num_vertices();
    const num_faces = mesh.num_faces();

    // -- XML header and VTKFile element --
    try writer.writeAll("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    try writer.writeAll("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    try writer.writeAll("  <UnstructuredGrid>\n");
    try writer.print("    <Piece NumberOfPoints=\"{d}\" NumberOfCells=\"{d}\">\n", .{ num_vertices, num_faces });

    // -- PointData (0-forms: one value per vertex) --
    if (point_data.len > 0) {
        try writer.writeAll("      <PointData>\n");
        for (point_data) |pd| {
            std.debug.assert(pd.values.len == num_vertices);
            try writeDataArray(writer, pd.name, pd.values);
        }
        try writer.writeAll("      </PointData>\n");
    }

    // -- CellData (1-forms or 2-forms mapped to faces) --
    if (cell_data.len > 0) {
        try writer.writeAll("      <CellData>\n");
        for (cell_data) |cd| {
            std.debug.assert(cd.values.len == num_faces);
            try writeDataArray(writer, cd.name, cd.values);
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

    // Connectivity: vertex indices for each triangle
    try writer.writeAll("        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");
    {
        const face_verts = mesh.simplices(2).items(.vertices);
        for (face_verts) |verts| {
            try writer.print("          {d} {d} {d}\n", .{ verts[0], verts[1], verts[2] });
        }
    }
    try writer.writeAll("        </DataArray>\n");

    // Offsets: cumulative vertex count per cell (3, 6, 9, ...)
    try writer.writeAll("        <DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
    {
        try writer.writeAll("          ");
        for (0..num_faces) |f| {
            if (f > 0) try writer.writeByte(' ');
            try writer.print("{d}", .{(f + 1) * 3});
        }
        try writer.writeByte('\n');
    }
    try writer.writeAll("        </DataArray>\n");

    // Types: all triangles (VTK type 5)
    try writer.writeAll("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    {
        try writer.writeAll("          ");
        for (0..num_faces) |f| {
            if (f > 0) try writer.writeByte(' ');
            try writer.print("{d}", .{vtk_triangle});
        }
        try writer.writeByte('\n');
    }
    try writer.writeAll("        </DataArray>\n");

    try writer.writeAll("      </Cells>\n");

    // -- Close elements --
    try writer.writeAll("    </Piece>\n");
    try writer.writeAll("  </UnstructuredGrid>\n");
    try writer.writeAll("</VTKFile>\n");
}

/// Write a single `<DataArray>` element with Float64 scalar data.
fn writeDataArray(writer: anytype, name: []const u8, values: []const f64) !void {
    try writer.print("        <DataArray type=\"Float64\" Name=\"{s}\" format=\"ascii\">\n", .{name});
    try writer.writeAll("          ");
    for (values, 0..) |v, i| {
        if (i > 0) try writer.writeByte(' ');
        try writer.print("{e}", .{v});
    }
    try writer.writeByte('\n');
    try writer.writeAll("        </DataArray>\n");
}

// ═══════════════════════════════════════════════════════════════════════════
// Field projection: 1-forms (edges) → per-face averages
// ═══════════════════════════════════════════════════════════════════════════

/// Project a 1-form (one value per edge) to a per-face scalar by averaging
/// the absolute values of the three edge coefficients on each face.
///
/// VTK CellData requires one value per cell (face). A 1-form lives on edges,
/// so direct export is not possible. This projection averages |E| over each
/// triangle's three boundary edges, giving a per-face proxy for field intensity.
///
/// Uses the ∂₂ boundary matrix (face → edge incidence) to find each face's edges.
pub fn project_edges_to_faces(
    allocator: std.mem.Allocator,
    comptime embedding_dimension: usize,
    comptime topological_dimension: usize,
    mesh: topology.Mesh(embedding_dimension, topological_dimension),
    edge_values: []const f64,
) ![]f64 {
    const num_faces = mesh.num_faces();
    const num_edges = mesh.num_edges();
    std.debug.assert(edge_values.len == num_edges);

    const face_values = try allocator.alloc(f64, num_faces);

    const boundary_2 = mesh.boundary(2);
    for (0..num_faces) |f| {
        const start = boundary_2.row_ptr[f];
        const end = boundary_2.row_ptr[f + 1];
        const edge_indices = boundary_2.col_idx[start..end];

        var sum: f64 = 0;
        for (edge_indices) |edge_idx| {
            sum += @abs(edge_values[edge_idx]);
        }
        face_values[f] = sum / @as(f64, @floatFromInt(edge_indices.len));
    }

    return face_values;
}

/// Write a VTK snapshot of electromagnetic fields (E and B) from a Maxwell state.
///
/// E is a primal 1-form (per-edge) and B is a primal 2-form (per-face).
/// Since VTK CellData requires per-face values, E is projected onto faces
/// by averaging |E| over each triangle's boundary edges.
///
/// The output .vtu file contains two CellData arrays:
///   - "E_intensity": per-face average of |E| on boundary edges
///   - "B_flux": per-face magnetic flux (direct from B cochain)
pub fn write_fields(
    allocator: std.mem.Allocator,
    writer: anytype,
    comptime embedding_dimension: usize,
    comptime topological_dimension: usize,
    mesh: topology.Mesh(embedding_dimension, topological_dimension),
    e_values: []const f64,
    b_values: []const f64,
) !void {
    const e_projected = try project_edges_to_faces(allocator, embedding_dimension, topological_dimension, mesh, e_values);
    defer allocator.free(e_projected);

    const cell_data = [_]DataArraySlice{
        .{ .name = "E_intensity", .values = e_projected },
        .{ .name = "B_flux", .values = b_values },
    };
    try write(writer, embedding_dimension, topological_dimension, mesh, &.{}, &cell_data);
}

// ═══════════════════════════════════════════════════════════════════════════
// Time-series snapshots
// ═══════════════════════════════════════════════════════════════════════════

/// A single entry in a PVD time-series collection.
pub const PvdEntry = struct {
    /// Simulation time for this snapshot (used by ParaView's animation timeline).
    timestep: f64,
    /// Relative path to the .vtu file (e.g., "output_0000.vtu").
    filename: []const u8,
};

/// Maximum length of a snapshot filename produced by `snapshot_filename`.
/// Sized for base names up to 200 chars + "_" + 10-digit u32 + ".vtu" = 215.
/// Rounded up to 224 for alignment headroom.
pub const max_snapshot_filename_length = 224;

/// Format a snapshot filename: "{base}_{step:0>4}.vtu".
///
/// Returns a slice into `buf` containing the formatted name. The step index
/// is zero-padded to 4 digits (0000–9999). Steps beyond 9999 use more digits.
///
/// Example: `snapshot_filename(&buf, "output", 42)` → `"output_0042.vtu"`
pub fn snapshot_filename(
    buf: *[max_snapshot_filename_length]u8,
    base_name: []const u8,
    step: u32,
) []const u8 {
    std.debug.assert(base_name.len > 0);
    std.debug.assert(base_name.len <= 200);
    return std.fmt.bufPrint(buf, "{s}_{d:0>4}.vtu", .{ base_name, step }) catch unreachable;
}

/// Write a PVD (ParaView Data) collection file referencing a set of .vtu snapshots.
///
/// The PVD format is a lightweight XML wrapper that tells ParaView which .vtu
/// files belong to a time series and what simulation time each represents.
/// ParaView loads the collection as an animation, stepping through snapshots
/// in timestep order.
///
/// Reference: VTK File Formats specification, §19.3 "XML File Formats"
/// (Kitware, The VTK User's Guide, 11th Edition).
pub fn write_pvd(writer: anytype, entries: []const PvdEntry) !void {
    std.debug.assert(entries.len > 0);

    try writer.writeAll("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    try writer.writeAll("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    try writer.writeAll("  <Collection>\n");

    for (entries) |entry| {
        try writer.print("    <DataSet timestep=\"{e}\" file=\"{s}\"/>\n", .{
            entry.timestep,
            entry.filename,
        });
    }

    try writer.writeAll("  </Collection>\n");
    try writer.writeAll("</VTKFile>\n");
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

    try write(output.writer(allocator), 2, 2, mesh, &.{}, &.{});

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

    try write(output.writer(allocator), 2, 2, mesh, &.{}, &.{});

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
    try write(output.writer(allocator), 2, 2, mesh, &pd, &.{});

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
    try write(output.writer(allocator), 2, 2, mesh, &.{}, &cd);

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
    try write(output.writer(allocator), 2, 2, mesh, &pd, &cd);

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

    try write(output.writer(allocator), 2, 2, mesh, &.{}, &.{});

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

    try write(output.writer(allocator), 2, 2, mesh, &.{}, &.{});

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

// ═══════════════════════════════════════════════════════════════════════════
// Time-series tests
// ═══════════════════════════════════════════════════════════════════════════

test "snapshot_filename formats zero-padded index" {
    var buf: [max_snapshot_filename_length]u8 = undefined;

    const name0 = snapshot_filename(&buf, "output", 0);
    try testing.expectEqualStrings("output_0000.vtu", name0);

    const name42 = snapshot_filename(&buf, "output", 42);
    try testing.expectEqualStrings("output_0042.vtu", name42);

    const name9999 = snapshot_filename(&buf, "output", 9999);
    try testing.expectEqualStrings("output_9999.vtu", name9999);
}

test "snapshot_filename handles step beyond 4 digits" {
    var buf: [max_snapshot_filename_length]u8 = undefined;
    const name = snapshot_filename(&buf, "sim", 12345);
    try testing.expectEqualStrings("sim_12345.vtu", name);
}

test "pvd output contains valid XML structure" {
    const allocator = testing.allocator;
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const entries = [_]PvdEntry{
        .{ .timestep = 0.0, .filename = "field_0000.vtu" },
        .{ .timestep = 0.5, .filename = "field_0001.vtu" },
        .{ .timestep = 1.0, .filename = "field_0002.vtu" },
    };

    try write_pvd(output.writer(allocator), &entries);

    const xml = output.items;

    // XML declaration
    try testing.expect(std.mem.startsWith(u8, xml, "<?xml version=\"1.0\""));
    // PVD root element (Collection type, not UnstructuredGrid)
    try testing.expect(std.mem.indexOf(u8, xml, "<VTKFile type=\"Collection\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "<Collection>") != null);
    // All three DataSet entries present with correct filenames
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0001.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0002.vtu\"") != null);
    // Closing elements
    try testing.expect(std.mem.indexOf(u8, xml, "</Collection>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "</VTKFile>") != null);
}

test "pvd entries reference correct timesteps" {
    const allocator = testing.allocator;
    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    const dt = 0.025;
    const entries = [_]PvdEntry{
        .{ .timestep = 0.0 * dt, .filename = "wave_0000.vtu" },
        .{ .timestep = 1.0 * dt, .filename = "wave_0001.vtu" },
        .{ .timestep = 2.0 * dt, .filename = "wave_0002.vtu" },
    };

    try write_pvd(output.writer(allocator), &entries);

    const xml = output.items;

    // Each DataSet element is self-closing and references both timestep and file.
    // Verify the timestep values appear in the output (formatted as scientific notation).
    try testing.expect(std.mem.indexOf(u8, xml, "timestep=\"0e0\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"wave_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"wave_0002.vtu\"") != null);
}

// ═══════════════════════════════════════════════════════════════════════════
// Field projection and write_fields tests (#42)
// ═══════════════════════════════════════════════════════════════════════════

test "project_edges_to_faces averages absolute edge values per face" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // Set each edge value to its index for distinct, verifiable values.
    const edge_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(edge_values);
    for (edge_values, 0..) |*v, i| {
        v.* = @as(f64, @floatFromInt(i)) + 1.0;
    }

    const face_values = try project_edges_to_faces(allocator, 2, 2, mesh, edge_values);
    defer allocator.free(face_values);

    // Each face should have exactly 3 edges (triangles).
    // Verify: for each face, the projected value equals the mean of its 3 edge |values|.
    const boundary_2 = mesh.boundary(2);
    for (0..mesh.num_faces()) |f| {
        const start = boundary_2.row_ptr[f];
        const end = boundary_2.row_ptr[f + 1];
        const edge_indices = boundary_2.col_idx[start..end];
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
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 1, 1, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // All negative edge values — projection should use absolute values.
    const edge_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(edge_values);
    for (edge_values) |*v| v.* = -3.0;

    const face_values = try project_edges_to_faces(allocator, 2, 2, mesh, edge_values);
    defer allocator.free(face_values);

    for (face_values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 3.0), v, 1e-15);
    }
}

test "write_fields produces vtu with E_intensity and B_flux CellData" {
    const allocator = testing.allocator;
    var mesh = try topology.Mesh(2, 2).uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // Create E (1-form, per-edge) and B (2-form, per-face).
    const e_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(e_values);
    for (e_values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) * 0.1;

    const b_values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(b_values);
    for (b_values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) * 0.5;

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write_fields(allocator, output.writer(allocator), 2, 2, mesh, e_values, b_values);

    const xml = output.items;

    // Must contain both field arrays in CellData.
    try testing.expect(std.mem.indexOf(u8, xml, "<CellData>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"E_intensity\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "Name=\"B_flux\"") != null);

    // B_flux values should round-trip exactly.
    const parsed_b = try parseDataArray(allocator, xml, "B_flux");
    defer allocator.free(parsed_b);
    try testing.expectEqual(b_values.len, parsed_b.len);
    for (b_values, parsed_b) |expected, actual| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }

    // E_intensity should have one value per face (projected from edges).
    const parsed_e = try parseDataArray(allocator, xml, "E_intensity");
    defer allocator.free(parsed_e);
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(parsed_e.len)));
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
        try write(output.writer(allocator), 2, 2, mesh, &pd, &.{});

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
        try write(output.writer(allocator), 2, 2, mesh, &.{}, &cd);

        const parsed = try parseDataArray(allocator, output.items, "random_flux");
        defer allocator.free(parsed);

        try testing.expectEqual(values.len, parsed.len);
        for (values, parsed) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-15);
        }
    }
}

test "snapshot_filename and write_pvd compose for time-series workflow" {
    // Simulate a 3-step time series: generate filenames, then write PVD.
    const allocator = testing.allocator;
    const num_steps = 3;
    const dt = 0.1;

    var filenames: [num_steps][max_snapshot_filename_length]u8 = undefined;
    var entries: [num_steps]PvdEntry = undefined;

    for (0..num_steps) |i| {
        const name = snapshot_filename(&filenames[i], "field", @intCast(i));
        entries[i] = .{
            .timestep = @as(f64, @floatFromInt(i)) * dt,
            .filename = name,
        };
    }

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try write_pvd(output.writer(allocator), &entries);

    const xml = output.items;

    // Verify the generated filenames appear in the PVD output.
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0001.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0002.vtu\"") != null);
    // Verify it is well-formed XML
    try testing.expect(std.mem.indexOf(u8, xml, "<VTKFile type=\"Collection\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "</VTKFile>") != null);
}
