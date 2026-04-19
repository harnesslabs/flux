//! FEEC weak-form assembly surfaces for the Laplace-de Rham operator.
//!
//! The strong DEC operator `Δₖ` acts on cochain values. In FEEC, the natural
//! assembled sparse object is the weak bilinear form matrix `Sₖ`, related by
//!
//!   Sₖ ω = Mₖ (Δₖ ω)
//!
//! where `Mₖ` is the k-form mass matrix. This module owns that weak assembly
//! surface.

const std = @import("std");
const sparse = @import("../../math/sparse.zig");

pub fn supportsExplicitStiffness(comptime MeshType: type, comptime k: comptime_int) bool {
    return switch (k) {
        0 => true,
        1 => MeshType.topological_dimension == 3,
        else => false,
    };
}

pub fn assemble_stiffness(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const MeshType = @TypeOf(mesh.*);
    const n = MeshType.topological_dimension;

    comptime {
        if (k < 0 or k > n) {
            @compileError(std.fmt.comptimePrint(
                "no Laplacian degree {d} on a {d}-dimensional mesh",
                .{ k, n },
            ));
        }
    }

    return switch (k) {
        0 => assemble_zero_form_stiffness(allocator, mesh),
        1 => switch (n) {
            3 => assemble_one_form_stiffness(allocator, mesh),
            else => @compileError("explicit FEEC stiffness assembly currently supports degree 1 only on tetrahedral meshes"),
        },
        else => @compileError("explicit FEEC stiffness assembly is not implemented for this degree"),
    };
}

fn assemble_zero_form_stiffness(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    return assemble_d_term_stiffness(0, allocator, mesh);
}

fn assemble_one_form_stiffness(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    var curl_curl = try assemble_d_term_stiffness(1, allocator, mesh);
    errdefer curl_curl.deinit(allocator);

    var grad_div = try assemble_one_form_grad_div_stiffness(allocator, mesh);
    errdefer grad_div.deinit(allocator);

    var assembler = sparse.TripletAssembler(f64).init(mesh.num_edges(), mesh.num_edges());
    defer assembler.deinit(allocator);

    try appendMatrixTriplets(&assembler, allocator, curl_curl);
    try appendMatrixTriplets(&assembler, allocator, grad_div);

    curl_curl.deinit(allocator);
    grad_div.deinit(allocator);
    return assembler.build(allocator);
}

fn assemble_d_term_stiffness(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const derivative = mesh.boundary(k + 1);
    const mass = mesh.whitney_mass(k + 1);
    const dof_count = mesh.num_cells(k);

    var assembler = sparse.TripletAssembler(f64).init(dof_count, dof_count);
    defer assembler.deinit(allocator);

    for (0..mass.n_rows) |cell_i| {
        const incidence_i = derivative.row(@intCast(cell_i));
        const mass_row = mass.row(@intCast(cell_i));

        for (mass_row.cols, mass_row.vals) |cell_j, mass_ij| {
            const incidence_j = derivative.row(cell_j);

            for (incidence_i.cols, 0..) |vertex_i, entry_i| {
                const sign_i = incidence_i.sign(entry_i);
                const left = @as(f64, @floatFromInt(sign_i)) * mass_ij;
                for (incidence_j.cols, 0..) |vertex_j, entry_j| {
                    const sign_j = incidence_j.sign(entry_j);
                    const contribution = left * @as(f64, @floatFromInt(sign_j));
                    try assembler.addEntry(allocator, vertex_i, vertex_j, contribution);
                }
            }
        }
    }

    return assembler.build(allocator);
}

fn assemble_one_form_grad_div_stiffness(
    allocator: std.mem.Allocator,
    mesh: anytype,
) !sparse.CsrMatrix(f64) {
    const edge_boundary = mesh.boundary(1);
    const mass = mesh.whitney_mass(1);
    const dual_volumes = mesh.vertices.slice().items(.dual_volume);
    const edge_count = mesh.num_edges();
    const vertex_count = mesh.num_vertices();

    const incidence_counts = try allocator.alloc(u32, vertex_count);
    defer allocator.free(incidence_counts);
    @memset(incidence_counts, 0);

    for (0..edge_count) |edge_idx| {
        const row = edge_boundary.row(@intCast(edge_idx));
        for (row.cols) |vertex_idx| {
            incidence_counts[vertex_idx] += 1;
        }
    }

    const vertex_offsets = try allocator.alloc(u32, @as(usize, vertex_count) + 1);
    defer allocator.free(vertex_offsets);
    vertex_offsets[0] = 0;
    for (0..vertex_count) |vertex_idx| {
        vertex_offsets[vertex_idx + 1] = vertex_offsets[vertex_idx] + incidence_counts[vertex_idx];
    }

    const total_incidence = vertex_offsets[vertex_count];
    const incident_edges = try allocator.alloc(u32, total_incidence);
    defer allocator.free(incident_edges);
    const incident_signs = try allocator.alloc(i8, total_incidence);
    defer allocator.free(incident_signs);

    const fill_offsets = try allocator.dupe(u32, vertex_offsets[0..vertex_count]);
    defer allocator.free(fill_offsets);

    for (0..edge_count) |edge_idx| {
        const row = edge_boundary.row(@intCast(edge_idx));
        for (row.cols, 0..) |vertex_idx, local_idx| {
            const write_idx = fill_offsets[vertex_idx];
            incident_edges[write_idx] = @intCast(edge_idx);
            incident_signs[write_idx] = row.sign(local_idx);
            fill_offsets[vertex_idx] += 1;
        }
    }

    var assembler = sparse.TripletAssembler(f64).init(edge_count, edge_count);
    defer assembler.deinit(allocator);

    var edge_weights = std.AutoHashMap(u32, f64).init(allocator);
    defer edge_weights.deinit();
    var nonzero_entries = std.ArrayListUnmanaged(struct {
        edge_idx: u32,
        weight: f64,
    }){};
    defer nonzero_entries.deinit(allocator);

    for (0..vertex_count) |vertex_idx| {
        edge_weights.clearRetainingCapacity();
        nonzero_entries.clearRetainingCapacity();

        const inv_dual_volume = 1.0 / dual_volumes[vertex_idx];
        std.debug.assert(dual_volumes[vertex_idx] != 0.0);

        const start = vertex_offsets[vertex_idx];
        const end = vertex_offsets[vertex_idx + 1];
        for (start..end) |entry_idx| {
            const incident_edge = incident_edges[entry_idx];
            const sign: f64 = @floatFromInt(incident_signs[entry_idx]);
            const mass_row = mass.row(incident_edge);
            for (mass_row.cols, mass_row.vals) |edge_j, mass_value| {
                const gop = try edge_weights.getOrPut(edge_j);
                if (!gop.found_existing) {
                    gop.value_ptr.* = 0.0;
                }
                gop.value_ptr.* += sign * mass_value;
            }
        }

        var iterator = edge_weights.iterator();
        while (iterator.next()) |entry| {
            try nonzero_entries.append(allocator, .{
                .edge_idx = entry.key_ptr.*,
                .weight = entry.value_ptr.*,
            });
        }

        for (nonzero_entries.items) |left| {
            const scaled_left = inv_dual_volume * left.weight;
            for (nonzero_entries.items) |right| {
                try assembler.addEntry(
                    allocator,
                    left.edge_idx,
                    right.edge_idx,
                    scaled_left * right.weight,
                );
            }
        }
    }

    return assembler.build(allocator);
}

fn appendMatrixTriplets(
    assembler: *sparse.TripletAssembler(f64),
    allocator: std.mem.Allocator,
    matrix: sparse.CsrMatrix(f64),
) !void {
    for (0..matrix.n_rows) |row_idx| {
        const row = matrix.row(@intCast(row_idx));
        for (row.cols, row.vals) |col_idx, value| {
            try assembler.addEntry(allocator, @intCast(row_idx), col_idx, value);
        }
    }
}
