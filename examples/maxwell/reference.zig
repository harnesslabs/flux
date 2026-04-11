const std = @import("std");
const flux = @import("flux");
const common = @import("examples_common");
const runtime = @import("runtime.zig");

const cochain = flux.forms;
const operator_context_mod = flux.operators.context;

pub fn project_te10_e(mesh: *const runtime.Mesh2D, values: []f64, t: f64, domain_length: f64) void {
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    const k = std.math.pi / domain_length;
    const omega = k;
    for (values, edge_verts) |*val, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const dy = p1[1] - p0[1];
        val.* = @sin(k * mx) * @sin(omega * t) * dy;
    }
}

pub fn project_te10_b(mesh: *const runtime.Mesh2D, values: []f64, t: f64, domain_length: f64) void {
    const simplex_2 = mesh.simplices(2);
    const face_verts = simplex_2.items(.vertices);
    const face_areas = simplex_2.items(.volume);
    const coords = mesh.vertices.slice().items(.coords);
    const k = std.math.pi / domain_length;
    const omega = k;
    for (values, face_verts, face_areas) |*val, verts, area| {
        const cx = (coords[verts[0]][0] + coords[verts[1]][0] + coords[verts[2]][0]) / 3.0;
        val.* = @cos(k * cx) * @cos(omega * t) * area;
    }
}

fn discrete_energy_from_arrays(allocator: std.mem.Allocator, mesh: *const runtime.Mesh2D, e_vals: []const f64, b_vals: []const f64) !f64 {
    const operator_context = try operator_context_mod.OperatorContext(runtime.Mesh2D).init(allocator, mesh);
    defer operator_context.deinit();
    _ = try operator_context.hodgeStar(1);
    _ = try operator_context.hodgeStar(2);

    var E = try runtime.MaxwellState2D.OneForm.init(allocator, mesh);
    defer E.deinit(allocator);
    @memcpy(E.values, e_vals);

    var B = try runtime.MaxwellState2D.TwoForm.init(allocator, mesh);
    defer B.deinit(allocator);
    @memcpy(B.values, b_vals);

    const state = struct {
        E: runtime.MaxwellState2D.OneForm,
        B: runtime.MaxwellState2D.TwoForm,
        operators: *operator_context_mod.OperatorContext(runtime.Mesh2D),
    }{ .E = E, .B = B, .operators = operator_context };
    return runtime.electromagnetic_energy(allocator, state);
}

pub fn run_cavity_energy_drift(allocator: std.mem.Allocator, grid_n: u32, final_time: f64) !f64 {
    const domain_length: f64 = 1.0;
    const dt = 0.1 * (domain_length / @as(f64, @floatFromInt(grid_n)));
    const num_steps: u32 = @intFromFloat(@round(final_time / dt));

    var mesh = try runtime.Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);

    var state = try runtime.MaxwellState2D.init(allocator, &mesh);
    defer state.deinit(allocator);

    project_te10_b(&mesh, state.B.values, -dt / 2.0, domain_length);
    const energy_initial = try discrete_energy_from_arrays(allocator, &mesh, state.E.values, state.B.values);

    for (0..num_steps) |_| {
        try runtime.leapfrog_step(allocator, &state, dt);
        runtime.apply_pec_boundary(&state);
    }

    const energy_final = try discrete_energy_from_arrays(allocator, &mesh, state.E.values, state.B.values);
    return @abs(energy_final - energy_initial) / energy_initial;
}

pub fn compute_te10_eigenvalue(allocator: std.mem.Allocator, grid_n: u32, domain_length: f64) !f64 {
    var mesh = try runtime.Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);
    const operator_context = try operator_context_mod.OperatorContext(runtime.Mesh2D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.exteriorDerivative(cochain.Primal, 1);
    _ = try operator_context.hodgeStar(1);
    _ = try operator_context.hodgeStar(2);

    var E = try runtime.MaxwellState2D.OneForm.init(allocator, &mesh);
    defer E.deinit(allocator);
    project_te10_e(&mesh, E.values, domain_length / 2.0, domain_length);

    var dE = try (try operator_context.exteriorDerivative(cochain.Primal, 1)).apply(allocator, E);
    defer dE.deinit(allocator);
    var star_dE = try (try operator_context.hodgeStar(2)).apply(allocator, dE);
    defer star_dE.deinit(allocator);
    var numerator: f64 = 0.0;
    for (dE.values, star_dE.values) |de, sde| numerator += de * sde;

    var star_E = try (try operator_context.hodgeStar(1)).apply(allocator, E);
    defer star_E.deinit(allocator);
    var denominator: f64 = 0.0;
    for (E.values, star_E.values) |e, se| denominator += e * se;
    return numerator / denominator;
}

pub fn tm110AngularFrequency(width: f64, height: f64) f64 {
    const kx = std.math.pi / width;
    const ky = std.math.pi / height;
    return std.math.sqrt(kx * kx + ky * ky);
}

pub fn project_tm110_e(mesh: *const runtime.Mesh3D, values: []f64, t: f64, width: f64, height: f64) void {
    const kx = std.math.pi / width;
    const ky = std.math.pi / height;
    const omega = tm110AngularFrequency(width, height);
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const my = 0.5 * (p0[1] + p1[1]);
        const dz = p1[2] - p0[2];
        value.* = @sin(kx * mx) * @sin(ky * my) * @sin(omega * t) * dz;
    }
}

pub fn project_tm110_potential(mesh: *const runtime.Mesh3D, values: []f64, t: f64, width: f64, height: f64) void {
    const kx = std.math.pi / width;
    const ky = std.math.pi / height;
    const omega = tm110AngularFrequency(width, height);
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const my = 0.5 * (p0[1] + p1[1]);
        const dz = p1[2] - p0[2];
        value.* = (1.0 / omega) * @sin(kx * mx) * @sin(ky * my) * @cos(omega * t) * dz;
    }
}

pub fn project_tm110_b(mesh: *const runtime.Mesh3D, values: []f64, t: f64, width: f64, height: f64) void {
    const kx = std.math.pi / width;
    const ky = std.math.pi / height;
    const omega = tm110AngularFrequency(width, height);
    const face_verts = mesh.simplices(2).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, face_verts) |*value, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const p2 = coords[verts[2]];
        const cx = (p0[0] + p1[0] + p2[0]) / 3.0;
        const cy = (p0[1] + p1[1] + p2[1]) / 3.0;
        const bx = (ky / omega) * @sin(kx * cx) * @cos(ky * cy) * @cos(omega * t);
        const by = -(kx / omega) * @cos(kx * cx) * @sin(ky * cy) * @cos(omega * t);
        const ea = [3]f64{ p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
        const eb = [3]f64{ p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
        const area_vector = [3]f64{
            0.5 * (ea[1] * eb[2] - ea[2] * eb[1]),
            0.5 * (ea[2] * eb[0] - ea[0] * eb[2]),
            0.5 * (ea[0] * eb[1] - ea[1] * eb[0]),
        };
        value.* = bx * area_vector[0] + by * area_vector[1];
    }
}

pub fn seedTm110Mode(allocator: std.mem.Allocator, state: anytype, dt: f64, width: f64, height: f64) !void {
    @memset(state.E.values, 0.0);
    var potential = try @TypeOf(state.*).OneForm.init(allocator, state.mesh);
    defer potential.deinit(allocator);
    project_tm110_potential(state.mesh, potential.values, -dt / 2.0, width, height);
    var exact_flux = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, potential);
    defer exact_flux.deinit(allocator);
    @memcpy(state.B.values, exact_flux.values);
}

pub fn writeSnapshot(allocator: std.mem.Allocator, writer: anytype, state: anytype) !void {
    try common.viz.writeProjectedTetFields(
        2,
        allocator,
        writer,
        state.mesh,
        .{
            .{ .name = "E_intensity", .kind = .edge_abs_mean, .values = state.E.values },
            .{ .name = "B_flux", .kind = .face_abs_mean, .values = state.B.values },
        },
    );
}

pub fn divergenceNorm3D(allocator: std.mem.Allocator, state: anytype) !f64 {
    var divergence = try (try state.operators.exteriorDerivative(cochain.Primal, 2)).apply(allocator, state.B);
    defer divergence.deinit(allocator);
    return std.math.sqrt(divergence.norm_squared());
}

pub fn compute_tm110_eigenvalue(allocator: std.mem.Allocator, nx: u32, ny: u32, nz: u32, width: f64, height: f64, depth: f64) !f64 {
    var mesh = try runtime.Mesh3D.uniform_tetrahedral_grid(allocator, nx, ny, nz, width, height, depth);
    defer mesh.deinit(allocator);
    const operator_context = try operator_context_mod.OperatorContext(runtime.Mesh3D).init(allocator, &mesh);
    defer operator_context.deinit();
    _ = try operator_context.exteriorDerivative(cochain.Primal, 1);
    _ = try operator_context.hodgeStar(1);
    _ = try operator_context.hodgeStar(2);

    var E = try runtime.MaxwellState3D.OneForm.init(allocator, &mesh);
    defer E.deinit(allocator);
    const omega = tm110AngularFrequency(width, height);
    project_tm110_e(&mesh, E.values, std.math.pi / (2.0 * omega), width, height);

    var dE = try (try operator_context.exteriorDerivative(cochain.Primal, 1)).apply(allocator, E);
    defer dE.deinit(allocator);
    var star_dE = try (try operator_context.hodgeStar(2)).apply(allocator, dE);
    defer star_dE.deinit(allocator);
    var numerator: f64 = 0.0;
    for (dE.values, star_dE.values) |lhs, rhs| numerator += lhs * rhs;

    var star_E = try (try operator_context.hodgeStar(1)).apply(allocator, E);
    defer star_E.deinit(allocator);
    var denominator: f64 = 0.0;
    for (E.values, star_E.values) |lhs, rhs| denominator += lhs * rhs;
    return numerator / denominator;
}
