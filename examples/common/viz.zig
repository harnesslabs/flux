const std = @import("std");

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
