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

/// Apply perfect electric conductor (PEC) boundary conditions.
///
/// Zeroes E on all mesh boundary edges. In the DEC formulation, E is a
/// primal 1-form (circulation along edges). PEC enforces that the tangential
/// electric field vanishes on the conducting boundary — equivalently, the
/// circulation of E along every boundary edge is zero.
///
/// Must be called after each Ampere step to enforce the constraint.
pub fn apply_pec_boundary(state: anytype) void {
    for (state.mesh.boundary_edges) |edge_idx| {
        state.E.values[edge_idx] = 0.0;
    }
}

/// Leapfrog integrator: advance (E, B) by one full timestep dt.
///
/// Uses the standard Yee-style staggered leapfrog where B lives at
/// half-integer times and E at integer times:
///
///   B^{n+1/2}  = B^{n-1/2} − dt · dE^n           (Faraday)
///   E^{n+1}    = E^n + dt · (★₁⁻¹ d ★₂ B^{n+1/2} − J)  (Ampere-Maxwell)
///
/// The caller stores B at the current half-step (n−1/2) and E at the
/// current integer step (n). After the call, B holds B^{n+1/2} and E
/// holds E^{n+1}. This preserves dB = 0 exactly at every half-step
/// because the Faraday update only adds a dd-exact term.
///
/// For source-free fields (J = 0), symplecticity of the leapfrog scheme
/// guarantees bounded energy oscillation — no secular drift.
pub fn leapfrog_step(
    allocator: std.mem.Allocator,
    state: anytype,
    dt: f64,
) !void {
    // Step 1: Faraday — advance B by a full dt at the half-step level.
    try faraday_step(allocator, state, dt);

    // Step 2: Ampere-Maxwell — advance E using the updated B.
    try ampere_step(allocator, state, dt);
}

/// Point dipole current source.
///
/// Models a localized sinusoidal current density concentrated on a single
/// edge — the edge whose midpoint is nearest to the given position. The
/// source produces J(t) = amplitude · sin(2π · frequency · t) on that
/// edge, with all other edges zero.
///
/// In 2D DEC, J is a primal 1-form with units of A/m integrated over the
/// edge length. The amplitude is divided by the edge length so that the
/// cochain value represents the correctly scaled line integral: the
/// circulation of J along the edge is amplitude · sin(...) / edge_length,
/// giving a density-like quantity consistent with the Ampere-Maxwell update.
pub fn PointDipole(comptime MeshType: type) type {
    return struct {
        const Self = @This();

        /// Index of the edge carrying the source.
        edge_index: u32,
        /// Length of the source edge (cached for scaling).
        edge_length: f64,
        /// Oscillation frequency in Hz.
        frequency: f64,
        /// Peak amplitude in amperes.
        amplitude: f64,

        /// Create a point dipole at the given position.
        ///
        /// Finds the edge whose midpoint is nearest to `position` and
        /// caches its index and length for repeated evaluation.
        pub fn init(mesh: *const MeshType, frequency: f64, amplitude: f64, position: [MeshType.dimension]f64) Self {
            const edge_slice = mesh.edges.slice();
            const edge_verts = edge_slice.items(.vertices);
            const lengths = edge_slice.items(.length);
            const coords = mesh.vertices.slice().items(.coords);

            var best_edge: u32 = 0;
            var best_distance_squared: f64 = std.math.inf(f64);

            for (0..mesh.num_edges()) |e| {
                const v0 = coords[edge_verts[e][0]];
                const v1 = coords[edge_verts[e][1]];

                // Edge midpoint.
                var midpoint: [MeshType.dimension]f64 = undefined;
                inline for (0..MeshType.dimension) |d| {
                    midpoint[d] = 0.5 * (v0[d] + v1[d]);
                }

                var dist_sq: f64 = 0;
                inline for (0..MeshType.dimension) |d| {
                    const diff = midpoint[d] - position[d];
                    dist_sq += diff * diff;
                }

                if (dist_sq < best_distance_squared) {
                    best_distance_squared = dist_sq;
                    best_edge = @intCast(e);
                }
            }

            return .{
                .edge_index = best_edge,
                .edge_length = lengths[best_edge],
                .frequency = frequency,
                .amplitude = amplitude,
            };
        }

        /// Evaluate the source at time `t` and write it into `J`.
        ///
        /// Sets J to zero everywhere except at the source edge, where
        /// J = (amplitude / edge_length) · sin(2π · frequency · t).
        pub fn apply(self: Self, J: anytype, t: f64) void {
            @memset(J.values, 0.0);
            J.values[self.edge_index] = (self.amplitude / self.edge_length) *
                @sin(2.0 * std.math.pi * self.frequency * t);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2);
const MaxwellState = State(Mesh2D);

// ── #37: Simulation state — timestep tracking ─────────────────────────

test "MaxwellState initializes with timestep zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try testing.expectEqual(@as(u64, 0), state.timestep);
}

test "leapfrog_step advances timestep" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try leapfrog_step(allocator, &state, 0.01);
    try testing.expectEqual(@as(u64, 1), state.timestep);

    try leapfrog_step(allocator, &state, 0.01);
    try testing.expectEqual(@as(u64, 2), state.timestep);
}

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

// ── #35: Leapfrog integrator ────────────────────────────────────────────

test "leapfrog_step composes Faraday then Ampere" {
    // A single leapfrog step should equal a Faraday step followed by an
    // Ampere step with the same dt.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    // State A: will use leapfrog_step.
    var state_a = try MaxwellState.init(allocator, &mesh);
    defer state_a.deinit(allocator);

    // State B: will use manual Faraday + Ampere.
    var state_b = try MaxwellState.init(allocator, &mesh);
    defer state_b.deinit(allocator);

    // Same random initial conditions.
    var rng = std.Random.DefaultPrng.init(0xE0_1F_01);
    for (state_a.E.values, state_b.E.values) |*a, *b| {
        const val = rng.random().float(f64) * 2.0 - 1.0;
        a.* = val;
        b.* = val;
    }
    for (state_a.B.values, state_b.B.values) |*a, *b| {
        const val = rng.random().float(f64) * 2.0 - 1.0;
        a.* = val;
        b.* = val;
    }

    const dt: f64 = 0.01;

    // Leapfrog.
    try leapfrog_step(allocator, &state_a, dt);

    // Manual: Faraday then Ampere.
    try faraday_step(allocator, &state_b, dt);
    try ampere_step(allocator, &state_b, dt);

    for (state_a.E.values, state_b.E.values) |a, b| {
        try testing.expectApproxEqAbs(b, a, 1e-15);
    }
    for (state_a.B.values, state_b.B.values) |a, b| {
        try testing.expectApproxEqAbs(b, a, 1e-15);
    }
}

test "leapfrog_step on zero fields is a no-op" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try leapfrog_step(allocator, &state, 0.01);

    for (state.E.values) |v| try testing.expectEqual(@as(f64, 0), v);
    for (state.B.values) |v| try testing.expectEqual(@as(f64, 0), v);
}

test "leapfrog energy bounded over 100 source-free steps" {
    // For source-free Maxwell (J = 0), the leapfrog integrator is
    // symplectic: total discrete energy should oscillate but not grow
    // secularly. We verify that energy at every step stays within a
    // bounded factor of the initial energy.
    const allocator = testing.allocator;

    // Use a fine enough grid so that the CFL condition is satisfied
    // for the chosen dt.
    var mesh = try Mesh2D.uniform_grid(allocator, 8, 8, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Seed with random E; B starts at zero (consistent with dB = 0).
    // Using small amplitudes to stay well within the linear/stable regime.
    var rng = std.Random.DefaultPrng.init(0xE0_1F_02);
    for (state.E.values) |*v| v.* = (rng.random().float(f64) * 2.0 - 1.0) * 0.01;

    // Compute initial energy: ‖E‖² + ‖B‖².
    const energy_initial = state.E.norm_squared() + state.B.norm_squared();

    // dt chosen small relative to the mesh spacing for stability.
    const dt: f64 = 0.001;
    const num_steps: usize = 100;

    // Track max energy ratio seen.
    var energy_ratio_max: f64 = 1.0;

    for (0..num_steps) |_| {
        try leapfrog_step(allocator, &state, dt);

        const energy = state.E.norm_squared() + state.B.norm_squared();
        const ratio = energy / energy_initial;
        energy_ratio_max = @max(energy_ratio_max, ratio);
    }

    // Symplectic integrator: energy should not drift secularly.
    // Allow up to 10% oscillation — the leapfrog scheme conserves a
    // shadow Hamiltonian, so the true energy oscillates by O(dt²).
    try testing.expect(energy_ratio_max < 1.1);
}

// Note: the real dB = 0 structural invariant test is #39. On a 2D mesh,
// d₂ maps into a zero-dimensional space (no 3-cells), so dB = 0 is
// vacuously true and cannot be tested meaningfully here. The energy
// boundedness test above already covers leapfrog stability.

// ── #36: PEC boundary conditions ─────────────────────────────────────────

test "apply_pec_boundary zeroes E on all boundary edges" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Fill E with non-zero values.
    for (state.E.values, 0..) |*v, i| v.* = @as(f64, @floatFromInt(i)) + 1.0;

    apply_pec_boundary(&state);

    // Boundary edges must be zero.
    for (mesh.boundary_edges) |edge_idx| {
        try testing.expectEqual(@as(f64, 0), state.E.values[edge_idx]);
    }

    // At least some interior edges should remain non-zero.
    var interior_nonzero: usize = 0;
    for (state.E.values) |v| {
        if (v != 0.0) interior_nonzero += 1;
    }
    try testing.expect(interior_nonzero > 0);
}

test "apply_pec_boundary is idempotent" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    var rng = std.Random.DefaultPrng.init(0xBEC_01);
    for (state.E.values) |*v| v.* = rng.random().float(f64);

    apply_pec_boundary(&state);

    // Save state after first application.
    const snapshot = try allocator.alloc(f64, state.E.values.len);
    defer allocator.free(snapshot);
    @memcpy(snapshot, state.E.values);

    // Apply again — should be identical.
    apply_pec_boundary(&state);
    for (state.E.values, snapshot) |actual, expected| {
        try testing.expectEqual(expected, actual);
    }
}

test "boundary edge count matches expected for rectangular grid" {
    // An nx × ny rectangular grid has 2*(nx + ny) boundary edges on the
    // rectangle sides, plus nx*ny diagonal edges are always interior.
    // Boundary edges: bottom nx + top nx + left ny + right ny = 2*(nx+ny).
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 5, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    try testing.expectEqual(@as(usize, 2 * (5 + 4)), mesh.boundary_edges.len);
}

// ── #38: Point dipole source ─────────────────────────────────────────────

const Dipole = PointDipole(Mesh2D);

test "PointDipole finds nearest edge to given position" {
    const allocator = testing.allocator;
    // 4×4 grid on [0,4]×[0,4], edge spacing = 1.
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 4.0, 4.0);
    defer mesh.deinit(allocator);

    // Place dipole at (0.5, 0.0) — midpoint of the first horizontal edge.
    const dipole = Dipole.init(&mesh, 1.0, 1.0, .{ 0.5, 0.0 });

    // The selected edge should have its midpoint very close to (0.5, 0).
    const verts = mesh.edges.slice().items(.vertices)[dipole.edge_index];
    const coords = mesh.vertices.slice().items(.coords);
    const mx = 0.5 * (coords[verts[0]][0] + coords[verts[1]][0]);
    const my = 0.5 * (coords[verts[0]][1] + coords[verts[1]][1]);

    try testing.expectApproxEqAbs(@as(f64, 0.5), mx, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.0), my, 1e-12);
}

test "PointDipole.apply produces sinusoidal source on one edge" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 4.0, 4.0);
    defer mesh.deinit(allocator);

    const frequency: f64 = 2.0;
    const amplitude: f64 = 3.0;
    const dipole = Dipole.init(&mesh, frequency, amplitude, .{ 2.0, 2.0 });

    var J = try MaxwellState.OneForm.init(allocator, &mesh);
    defer J.deinit(allocator);

    // At t = 0, sin(0) = 0 → J should be zero everywhere.
    dipole.apply(&J, 0.0);
    for (J.values) |v| try testing.expectEqual(@as(f64, 0), v);

    // At t = 1/(4·frequency), sin(π/2) = 1 → peak value on source edge.
    const t_quarter = 1.0 / (4.0 * frequency);
    dipole.apply(&J, t_quarter);

    const expected_peak = amplitude / dipole.edge_length;
    try testing.expectApproxEqAbs(expected_peak, J.values[dipole.edge_index], 1e-12);

    // All other edges should be zero.
    for (J.values, 0..) |v, i| {
        if (i != dipole.edge_index) {
            try testing.expectEqual(@as(f64, 0), v);
        }
    }
}

test "PointDipole.apply at half period is zero" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const dipole = Dipole.init(&mesh, 5.0, 1.0, .{ 0.5, 0.5 });
    var J = try MaxwellState.OneForm.init(allocator, &mesh);
    defer J.deinit(allocator);

    // At t = 1/(2·frequency), sin(π) ≈ 0.
    dipole.apply(&J, 1.0 / (2.0 * 5.0));
    for (J.values) |v| {
        try testing.expectApproxEqAbs(@as(f64, 0), v, 1e-12);
    }
}
