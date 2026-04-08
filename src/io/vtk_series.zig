const std = @import("std");
const testing = std.testing;

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
    try testing.expect(std.mem.startsWith(u8, xml, "<?xml version=\"1.0\""));
    try testing.expect(std.mem.indexOf(u8, xml, "<VTKFile type=\"Collection\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "<Collection>") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0001.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0002.vtu\"") != null);
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
    try testing.expect(std.mem.indexOf(u8, xml, "timestep=\"0e0\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"wave_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"wave_0002.vtu\"") != null);
}

test "snapshot_filename and write_pvd compose for time-series workflow" {
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
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0000.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0001.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "file=\"field_0002.vtu\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "<VTKFile type=\"Collection\"") != null);
    try testing.expect(std.mem.indexOf(u8, xml, "</VTKFile>") != null);
}
