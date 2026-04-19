//! Local-to-global assembly for FEEC weak-form operators.

const std = @import("std");
const sparse = @import("../math/sparse.zig");

pub fn assemble(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
    kernel: anytype,
) !sparse.CsrMatrix(f64) {
    const MeshType = @TypeOf(mesh.*);
    const n = MeshType.topological_dimension;

    comptime {
        if (k <= 0 or k >= n) {
            @compileError("weak-form assembly is only defined for interior degrees 0 < k < n");
        }
        if (n > 3) {
            @compileError("weak-form assembly is currently implemented for topological_dimension <= 3");
        }
        if (!@hasDecl(@TypeOf(kernel), "degree")) {
            @compileError("weak-form kernels must declare their primal degree");
        }
        if (@TypeOf(kernel).degree != k) {
            @compileError("weak-form kernel degree does not match assembly degree");
        }
        if (!@hasDecl(@TypeOf(kernel), "localMatrix")) {
            @compileError("weak-form kernels must provide localMatrix(top_simplex_index)");
        }
    }

    if (comptime @hasField(@TypeOf(kernel), "mesh")) {
        std.debug.assert(kernel.mesh == mesh);
    }

    const simplex_count = mesh.num_cells(k);
    const top_simplex_count = mesh.num_cells(n);
    const simplex_vertices = mesh.simplices(k).items(.vertices);
    const top_simplex_vertices = mesh.simplices(n).items(.vertices);
    const local_faces = localFaces(n, k);

    var simplex_index = std.AutoHashMap([k + 1]u32, u32).init(allocator);
    defer simplex_index.deinit();
    for (simplex_vertices, 0..) |vertices, global_idx| {
        try simplex_index.put(vertices, @intCast(global_idx));
    }

    var assembler = sparse.TripletAssembler(f64).init(simplex_count, simplex_count);
    defer assembler.deinit(allocator);

    for (0..top_simplex_count) |top_idx| {
        const local_matrix = kernel.localMatrix(top_idx);
        const top_vertices = top_simplex_vertices[top_idx];

        var global_indices: [local_faces.len]u32 = undefined;
        var orientation_signs: [local_faces.len]i8 = undefined;
        for (local_faces, 0..) |local_face, local_face_idx| {
            const oriented_vertices = liftLocalFaceVertices(k, n, top_vertices, local_face);
            const canonical_key = canonicalizeVertices(k + 1, oriented_vertices);
            const global_idx = simplex_index.get(canonical_key).?;
            global_indices[local_face_idx] = global_idx;
            orientation_signs[local_face_idx] = orientationSign(k + 1, oriented_vertices, simplex_vertices[global_idx]);
        }

        for (0..local_faces.len) |i| {
            const sign_i: f64 = @floatFromInt(orientation_signs[i]);
            for (0..local_faces.len) |j| {
                const sign_j: f64 = @floatFromInt(orientation_signs[j]);
                try assembler.addEntry(
                    allocator,
                    global_indices[i],
                    global_indices[j],
                    sign_i * sign_j * local_matrix[i][j],
                );
            }
        }
    }

    return assembler.build(allocator);
}

pub fn localFaceCount(comptime n: comptime_int, comptime k: comptime_int) comptime_int {
    return choose(n + 1, k + 1);
}

pub fn localFaces(comptime n: comptime_int, comptime k: comptime_int) [localFaceCount(n, k)][k + 1]u8 {
    return switch (n) {
        2 => switch (k) {
            1 => .{
                .{ 0, 1 },
                .{ 0, 2 },
                .{ 1, 2 },
            },
            else => @compileError("unsupported local face degree for 2-simplex"),
        },
        3 => switch (k) {
            1 => .{
                .{ 0, 1 },
                .{ 0, 2 },
                .{ 0, 3 },
                .{ 1, 2 },
                .{ 1, 3 },
                .{ 2, 3 },
            },
            2 => .{
                .{ 0, 1, 2 },
                .{ 0, 1, 3 },
                .{ 0, 2, 3 },
                .{ 1, 2, 3 },
            },
            else => @compileError("unsupported local face degree for 3-simplex"),
        },
        else => @compileError("local face enumeration is only implemented for topological_dimension <= 3"),
    };
}

fn choose(comptime n: comptime_int, comptime k: comptime_int) comptime_int {
    if (k < 0 or k > n) return 0;
    if (k == 0 or k == n) return 1;
    var numerator: comptime_int = 1;
    var denominator: comptime_int = 1;
    const k_small = if (k < n - k) k else n - k;
    inline for (0..k_small) |i| {
        numerator *= (n - i);
        denominator *= (i + 1);
    }
    return @divExact(numerator, denominator);
}

fn liftLocalFaceVertices(
    comptime k: comptime_int,
    comptime n: comptime_int,
    top_vertices: [n + 1]u32,
    local_face: [k + 1]u8,
) [k + 1]u32 {
    var result: [k + 1]u32 = undefined;
    for (local_face, 0..) |local_vertex_idx, write_idx| {
        result[write_idx] = top_vertices[local_vertex_idx];
    }
    return result;
}

fn canonicalizeVertices(comptime len: comptime_int, vertices: [len]u32) [len]u32 {
    var result = vertices;
    inline for (1..len) |i| {
        var j = i;
        while (j > 0 and result[j - 1] > result[j]) : (j -= 1) {
            std.mem.swap(u32, &result[j - 1], &result[j]);
        }
    }
    return result;
}

fn orientationSign(comptime len: comptime_int, oriented_vertices: [len]u32, global_vertices: [len]u32) i8 {
    var permutation: [len]u8 = undefined;
    var used = [_]bool{false} ** len;
    for (oriented_vertices, 0..) |vertex, oriented_idx| {
        var found = false;
        for (global_vertices, 0..) |global_vertex, global_idx| {
            if (used[global_idx]) continue;
            if (vertex != global_vertex) continue;
            permutation[oriented_idx] = @intCast(global_idx);
            used[global_idx] = true;
            found = true;
            break;
        }
        std.debug.assert(found);
    }

    var inversion_count: u32 = 0;
    inline for (0..len) |i| {
        inline for (i + 1..len) |j| {
            if (permutation[i] > permutation[j]) inversion_count += 1;
        }
    }
    return if (inversion_count % 2 == 0) 1 else -1;
}
