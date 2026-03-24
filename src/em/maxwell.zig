//! Maxwell's equations on a discrete exterior calculus mesh.
//!
//! The electromagnetic field is represented as:
//!   - E ∈ Ω¹ (primal 1-form on edges) — electric field circulation
//!   - B ∈ Ω² (primal 2-form on faces) — magnetic flux
//!   - J ∈ Ω¹ (primal 1-form on edges) — current density source
//!
//! The semi-discrete Maxwell system in DEC form:
//!   ∂B/∂t = −dE          (Faraday)
//!   ∂E/∂t = ★₁⁻¹ d ★₂ B − J   (Ampere-Maxwell)
//!
//! The remaining Maxwell equations are constraints, not evolution equations:
//!   d★E = ρ/ε₀           (Gauss's law — divergence of E)
//!   dB  = 0              (no magnetic monopoles)
//!   ∂ρ/∂t + d★J = 0     (charge conservation)
//!
//! The Faraday update preserves dB = 0 exactly because d(dE) = 0
//! algebraically — the divergence-free constraint is structural, not
//! numerical. Gauss's law and charge conservation are satisfied if
//! the initial conditions and source J are consistent; they are not
//! actively enforced per timestep.
//!
//! Update functions are standalone (not methods on `State`) to keep
//! integration logic separate from state representation, per the
//! pluggable time-integrator horizon.

const std = @import("std");
const testing = std.testing;
const cochain = @import("../forms/cochain.zig");
const topology = @import("../topology/mesh.zig");
const ext = @import("../operators/exterior_derivative.zig");
const hs = @import("../operators/hodge_star.zig");

/// Electromagnetic field state on a 2D simplicial mesh.
///
/// Holds the three physical fields (E, B, J) as typed cochains whose
/// degrees are enforced at compile time. The mesh is borrowed — the
/// caller must keep it alive for the lifetime of the state.
pub fn State(comptime MeshType: type) type {
    comptime {
        if (!@hasDecl(MeshType, "dimension")) {
            @compileError("Maxwell State requires a Mesh type with a 'dimension' declaration");
        }
    }

    return struct {
        const Self = @This();

        pub const OneForm = cochain.Cochain(MeshType, 1, cochain.Primal);
        pub const TwoForm = cochain.Cochain(MeshType, 2, cochain.Primal);

        /// Electric field — primal 1-form (circulation along edges).
        E: OneForm,
        /// Magnetic flux — primal 2-form (flux through faces).
        B: TwoForm,
        /// Current density source — primal 1-form (same space as E).
        J: OneForm,
        /// The mesh this state is defined on.
        mesh: *const MeshType,

        /// Allocate a zero-initialized Maxwell state on the given mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var E = try OneForm.init(allocator, mesh);
            errdefer E.deinit(allocator);

            var B = try TwoForm.init(allocator, mesh);
            errdefer B.deinit(allocator);

            var J = try OneForm.init(allocator, mesh);
            errdefer J.deinit(allocator);

            return .{ .E = E, .B = B, .J = J, .mesh = mesh };
        }

        /// Free all field storage.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.J.deinit(allocator);
            self.B.deinit(allocator);
            self.E.deinit(allocator);
        }
    };
}

/// Faraday update: B += −dt · d(E).
///
/// Applies the discrete exterior derivative to E (a primal 1-form)
/// to get dE (a primal 2-form), scales by dt, and subtracts from B.
///
/// This preserves d₂B = 0 exactly because:
///   d₂(B − dt · dE) = d₂B − dt · d₂(dE) = d₂B − 0 = d₂B
/// The dd = 0 identity is algebraic, not numerical.
pub fn faraday_step(
    allocator: std.mem.Allocator,
    state: anytype,
    dt: f64,
) !void {
    var dE = try ext.exterior_derivative(allocator, state.E);
    defer dE.deinit(allocator);

    dE.scale(dt);
    state.B.sub(dE);
}

/// Ampere-Maxwell update: E += dt · (★₁⁻¹ d ★₂ B − J).
///
/// Computes the discrete curl of B via ★₁⁻¹ ∘ d ∘ ★₂, subtracts the
/// current source J, scales by dt, and adds to E.
///
/// The operator chain types:
///   ★₂: primal 2-form → dual 0-form
///   d:  dual 0-form   → dual 1-form
///   ★₁⁻¹: dual 1-form → primal 1-form
///
/// Uses direct diagonal application (like the Hodge Laplacian) rather
/// than `hodge_star_inverse`, which panics on degenerate edges. On
/// uniform grids, diagonal edges have dual_length = 0, making ★₁
/// singular there. The physically correct treatment is ★₁⁻¹(0) = 0
/// at those edges — they carry no field information.
pub fn ampere_step(
    allocator: std.mem.Allocator,
    state: anytype,
    dt: f64,
) !void {
    // Step 1: ★₂(B) — primal 2-form → dual 0-form.
    var star_B = try hs.hodge_star(allocator, state.B);
    defer star_B.deinit(allocator);

    // Step 2: d(★₂B) — dual 0-form → dual 1-form.
    var d_star_B = try ext.exterior_derivative(allocator, star_B);
    defer d_star_B.deinit(allocator);

    // Step 3: ★₁⁻¹(d(★₂B)) — dual 1-form → primal 1-form.
    // Apply the ★₁ inverse diagonal manually to handle degenerate edges
    // (dual_length = 0) with the pseudo-inverse: 0/0 → 0.
    const edge_slice = state.mesh.edges.slice();
    const lengths = edge_slice.items(.length);
    const dual_lengths = edge_slice.items(.dual_length);

    for (state.E.values, d_star_B.values, lengths, dual_lengths, state.J.values) |*e, dsb, len, dual_len, j| {
        const curl_component = if (dual_len == 0.0)
            // Degenerate edge: ★₁ is singular here, pseudo-inverse gives 0.
            0.0
        else
            (len / dual_len) * dsb;

        e.* += dt * (curl_component - j);
    }

    // Comptime type assertion: verify the operator chain types are correct.
    // B is a primal 2-form. ★₂B is a dual 0-form. d(★₂B) is a dual 1-form.
    // ★₁⁻¹ maps dual 1-form → primal 1-form = E's type.
    comptime {
        const BType = @TypeOf(state.B);
        std.debug.assert(BType.degree == 2 and BType.duality == cochain.Primal);
        std.debug.assert(@TypeOf(state.E).degree == 1 and @TypeOf(state.E).duality == cochain.Primal);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);
const MaxwellState = State(Mesh2D);

// ── #32: Field assignments ─────────────────────────────────────────────

test "MaxwellState fields have correct degrees" {
    comptime {
        try testing.expectEqual(1, MaxwellState.OneForm.degree);
        try testing.expectEqual(2, MaxwellState.TwoForm.degree);
        try testing.expect(MaxwellState.OneForm.duality == cochain.Primal);
        try testing.expect(MaxwellState.TwoForm.duality == cochain.Primal);
    }
}

test "MaxwellState initializes E, B, J to zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.E.values.len)));
    try testing.expectEqual(mesh.num_faces(), @as(u32, @intCast(state.B.values.len)));
    try testing.expectEqual(mesh.num_edges(), @as(u32, @intCast(state.J.values.len)));

    for (state.E.values) |v| try testing.expectEqual(@as(f64, 0), v);
    for (state.B.values) |v| try testing.expectEqual(@as(f64, 0), v);
    for (state.J.values) |v| try testing.expectEqual(@as(f64, 0), v);
}

test "MaxwellState holds mesh reference" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 2, 2, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expect(state.mesh == &mesh);
    try testing.expect(state.E.mesh == &mesh);
    try testing.expect(state.B.mesh == &mesh);
}

// ── #33: Faraday update ────────────────────────────────────────────────

test "faraday_step on zero E leaves B unchanged" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set B to non-trivial values.
    for (state.B.values, 0..) |*v, i| v.* = @floatFromInt(i);

    // Save original B.
    const original_B = try allocator.alloc(f64, state.B.values.len);
    defer allocator.free(original_B);
    @memcpy(original_B, state.B.values);

    try faraday_step(allocator, &state, 0.01);

    // E = 0 ⟹ dE = 0 ⟹ B unchanged.
    for (state.B.values, original_B) |b, orig| {
        try testing.expectApproxEqAbs(orig, b, 1e-15);
    }
}

test "faraday_step produces non-trivial update for non-zero E" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set E to a non-constant function so dE ≠ 0.
    var rng = std.Random.DefaultPrng.init(0xE0_FAD_01);
    for (state.E.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    try faraday_step(allocator, &state, 0.01);

    // B should have at least some non-zero values.
    var max_abs: f64 = 0;
    for (state.B.values) |v| {
        max_abs = @max(max_abs, @abs(v));
    }
    try testing.expect(max_abs > 1e-10);
}

test "faraday_step is B -= dt * d(E)" {
    // Verify the update formula directly: B_new = B_old - dt * d(E).
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xE0_FAD_02);
    for (state.E.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;
    for (state.B.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    // Compute expected: B_old - dt * d(E)
    const dt: f64 = 0.025;
    var dE = try ext.exterior_derivative(allocator, state.E);
    defer dE.deinit(allocator);

    const expected_B = try allocator.alloc(f64, state.B.values.len);
    defer allocator.free(expected_B);
    for (expected_B, state.B.values, dE.values) |*exp, b, de| {
        exp.* = b - dt * de;
    }

    try faraday_step(allocator, &state, dt);

    for (state.B.values, expected_B) |actual, expected| {
        try testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

// ── #34: Ampere-Maxwell update ─────────────────────────────────────────

test "ampere_step on zero B and zero J leaves E unchanged" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set E to non-trivial values.
    for (state.E.values, 0..) |*v, i| v.* = @floatFromInt(i);

    const original_E = try allocator.alloc(f64, state.E.values.len);
    defer allocator.free(original_E);
    @memcpy(original_E, state.E.values);

    try ampere_step(allocator, &state, 0.01);

    for (state.E.values, original_E) |e, orig| {
        try testing.expectApproxEqAbs(orig, e, 1e-15);
    }
}

test "ampere_step produces non-trivial update for non-zero B" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set B to a non-constant function.
    var rng = std.Random.DefaultPrng.init(0xE0_A0B_01);
    for (state.B.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    try ampere_step(allocator, &state, 0.01);

    // E should have at least some non-zero values.
    var max_abs: f64 = 0;
    for (state.E.values) |v| {
        max_abs = @max(max_abs, @abs(v));
    }
    try testing.expect(max_abs > 1e-10);
}

test "ampere_step is E += dt * (★₁⁻¹ d ★₂ B − J)" {
    // Verify the update formula by manually computing each operator step
    // and comparing to the result of ampere_step.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xE0_A0B_02);
    for (state.E.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;
    for (state.B.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;
    for (state.J.values) |*v| v.* = rng.random().float(f64) * 0.1 - 0.05;

    // Compute expected manually: E_old + dt * (★₁⁻¹(d(★₂(B))) − J)
    const dt: f64 = 0.025;

    // ★₂(B)
    var star_B = try hs.hodge_star(allocator, state.B);
    defer star_B.deinit(allocator);

    // d(★₂B)
    var d_star_B = try ext.exterior_derivative(allocator, star_B);
    defer d_star_B.deinit(allocator);

    // ★₁⁻¹(d(★₂B)) with pseudo-inverse for degenerate edges
    const edge_slice = mesh.edges.slice();
    const lengths = edge_slice.items(.length);
    const dual_lengths = edge_slice.items(.dual_length);

    const expected_E = try allocator.alloc(f64, state.E.values.len);
    defer allocator.free(expected_E);
    for (expected_E, state.E.values, d_star_B.values, lengths, dual_lengths, state.J.values) |*exp, e, dsb, len, dual_len, j| {
        const curl_component = if (dual_len == 0.0) 0.0 else (len / dual_len) * dsb;
        exp.* = e + dt * (curl_component - j);
    }

    try ampere_step(allocator, &state, dt);

    for (state.E.values, expected_E) |actual, expected| {
        try testing.expectApproxEqAbs(expected, actual, 1e-14);
    }
}

test "ampere_step with J source reduces E compared to without" {
    // A positive current source J should reduce the E update compared
    // to the source-free case.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    // State without source
    var state_no_source = try MaxwellState.init(allocator, &mesh);
    defer state_no_source.deinit(allocator);

    // State with source
    var state_with_source = try MaxwellState.init(allocator, &mesh);
    defer state_with_source.deinit(allocator);

    // Same initial B in both.
    var rng = std.Random.DefaultPrng.init(0xE0_A0B_03);
    for (state_no_source.B.values, state_with_source.B.values) |*b1, *b2| {
        const val = rng.random().float(f64) * 2.0 - 1.0;
        b1.* = val;
        b2.* = val;
    }

    // Add a uniform positive source to the second state.
    for (state_with_source.J.values) |*v| v.* = 1.0;

    const dt: f64 = 0.01;
    try ampere_step(allocator, &state_no_source, dt);
    try ampere_step(allocator, &state_with_source, dt);

    // E_with_source = E_no_source − dt * J, so their difference should be −dt * J.
    for (state_with_source.E.values, state_no_source.E.values) |e_src, e_nosrc| {
        try testing.expectApproxEqAbs(e_nosrc - dt * 1.0, e_src, 1e-15);
    }
}
