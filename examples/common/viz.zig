const std = @import("std");
const flux = @import("flux");

pub const TetProjectionKind = enum {
    edge_abs_mean,
    face_abs_mean,
};

pub const TetProjectionField = struct {
    name: []const u8,
    kind: TetProjectionKind,
    values: []const f64,
};

pub fn projectEdgesToTets(
    allocator: std.mem.Allocator,
    mesh: anytype,
    edge_values: []const f64,
) ![]f64 {
    const output = try allocator.alloc(f64, mesh.num_tets());
    errdefer allocator.free(output);
    const touched = try allocator.alloc(bool, mesh.num_edges());
    defer allocator.free(touched);

    for (0..mesh.num_tets()) |tet_idx_usize| {
        @memset(touched, false);
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        var count: u32 = 0;

        for (tet_row.cols) |face_idx| {
            const face_row = mesh.boundary(2).row(face_idx);
            for (face_row.cols) |edge_idx| {
                if (touched[edge_idx]) continue;
                touched[edge_idx] = true;
                sum += @abs(edge_values[edge_idx]);
                count += 1;
            }
        }

        std.debug.assert(count > 0);
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(count));
    }

    return output;
}

pub fn projectFacesToTets(
    allocator: std.mem.Allocator,
    mesh: anytype,
    face_values: []const f64,
) ![]f64 {
    const output = try allocator.alloc(f64, mesh.num_tets());
    errdefer allocator.free(output);

    for (0..mesh.num_tets()) |tet_idx_usize| {
        const tet_row = mesh.boundary(3).row(@intCast(tet_idx_usize));
        var sum: f64 = 0.0;
        for (tet_row.cols) |face_idx| {
            sum += @abs(face_values[face_idx]);
        }
        output[tet_idx_usize] = sum / @as(f64, @floatFromInt(tet_row.cols.len));
    }

    return output;
}

pub fn writeProjectedTetFields(
    comptime field_count: usize,
    allocator: std.mem.Allocator,
    writer: anytype,
    mesh: anytype,
    fields: [field_count]TetProjectionField,
) !void {
    var projected_values: [field_count][]f64 = undefined;
    defer inline for (0..field_count) |field_idx| {
        allocator.free(projected_values[field_idx]);
    };

    var cell_data: [field_count]flux.io.DataArraySlice = undefined;
    inline for (fields, 0..) |field, field_idx| {
        projected_values[field_idx] = switch (field.kind) {
            .edge_abs_mean => try projectEdgesToTets(allocator, mesh, field.values),
            .face_abs_mean => try projectFacesToTets(allocator, mesh, field.values),
        };
        cell_data[field_idx] = .{
            .name = field.name,
            .values = projected_values[field_idx],
        };
    }

    try flux.io.write(writer, mesh.*, &.{}, &cell_data);
}

test "writeProjectedTetFields writes projected cell arrays" {
    const allocator = std.testing.allocator;
    const Mesh3D = flux.topology.Mesh(3, 3);

    var mesh = try Mesh3D.cartesian(allocator, .{ 1, 1, 1 }, .{ 1.0, 1.0, 1.0 });
    defer mesh.deinit(allocator);

    const edge_values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(edge_values);
    @memset(edge_values, 1.0);

    const face_values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(face_values);
    @memset(face_values, 2.0);

    var output = std.ArrayListUnmanaged(u8){};
    defer output.deinit(allocator);

    try writeProjectedTetFields(
        2,
        allocator,
        output.writer(allocator),
        &mesh,
        .{
            .{ .name = "edge_field", .kind = .edge_abs_mean, .values = edge_values },
            .{ .name = "face_field", .kind = .face_abs_mean, .values = face_values },
        },
    );

    try std.testing.expect(std.mem.indexOf(u8, output.items, "edge_field") != null);
    try std.testing.expect(std.mem.indexOf(u8, output.items, "face_field") != null);
}
