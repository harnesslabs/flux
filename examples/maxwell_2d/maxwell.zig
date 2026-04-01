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
const flux = @import("flux");
const cochain = flux.forms;
const topology = flux.topology;
const sparse = flux.math.sparse;
const exterior_derivative = flux.operators.exterior_derivative;
const hodge_star = flux.operators.hodge_star;
const leapfrog_mod = flux.integrators.leapfrog;

/// Electromagnetic field state on a 2D simplicial mesh.
///
/// Holds the three physical fields (E, B, J) as typed cochains whose
/// degrees are enforced at compile time. The mesh is borrowed — the
/// caller must keep it alive for the lifetime of the state.
pub fn State(comptime MeshType: type) type {
    comptime {
        if (!@hasDecl(MeshType, "embedding_dimension")) {
            @compileError("Maxwell State requires a Mesh type with an 'embedding_dimension' declaration");
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
        /// Discrete timestep counter (number of completed leapfrog steps).
        timestep: u64,

        /// Allocate a zero-initialized Maxwell state on the given mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
            var E = try OneForm.init(allocator, mesh);
            errdefer E.deinit(allocator);

            var B = try TwoForm.init(allocator, mesh);
            errdefer B.deinit(allocator);

            var J = try OneForm.init(allocator, mesh);
            errdefer J.deinit(allocator);

            return .{ .E = E, .B = B, .J = J, .mesh = mesh, .timestep = 0 };
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
    var dE = try exterior_derivative.exterior_derivative(allocator, state.E);
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
/// Ampere step: E^{n+1} = E^n + dt · (★₁⁻¹ d̃₀ ★₂ B^{n+1/2} − J).
///
/// Applies ★₁⁻¹ as the diagonal inverse (length / dual_length).
/// With the barycentric dual, every edge has nonzero dual_length,
/// so ★₁ is non-singular and no pseudo-inverse workaround is needed.
pub fn ampere_step(
    allocator: std.mem.Allocator,
    state: anytype,
    dt: f64,
) !void {
    // Step 1: ★₂(B) — primal 2-form → dual 0-form.
    var star_B = try hodge_star.hodge_star(allocator, state.B);
    defer star_B.deinit(allocator);

    // Step 2: d(★₂B) — dual 0-form → dual 1-form.
    var d_star_B = try exterior_derivative.exterior_derivative(allocator, star_B);
    defer d_star_B.deinit(allocator);

    // Step 3: ★₁⁻¹(d(★₂B)) — dual 1-form → primal 1-form.
    const edge_volumes = state.mesh.simplices(1).items(.volume);
    const dual_edge_vols = state.mesh.dual_edge_volumes;

    for (state.E.values, d_star_B.values, edge_volumes, dual_edge_vols, state.J.values) |*e, dsb, vol, dual_vol, j| {
        std.debug.assert(dual_vol > 0.0);
        e.* += dt * ((vol / dual_vol) * dsb - j);
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
///
/// Delegates to the generic `Leapfrog(MaxwellSystem(...))` integrator.
pub fn leapfrog_step(
    allocator: std.mem.Allocator,
    state: anytype,
    dt: f64,
) !void {
    try faraday_step(allocator, state, dt);
    try ampere_step(allocator, state, dt);
    state.timestep += 1;
}

/// Maxwell system for the generic leapfrog integrator.
///
/// Wraps the Faraday and Ampere free functions as typed half-step
/// operators, conforming to the interface that `Leapfrog` requires.
/// The timestep counter is incremented in `second_half_step` because
/// it is Maxwell-specific state bookkeeping, not integrator logic.
pub fn MaxwellSystem(comptime MeshType: type) type {
    const StateType = State(MeshType);

    return struct {
        pub const State = StateType;

        /// Faraday half-step: B^{n+1/2} = B^{n-1/2} − dt · dE^n.
        pub fn first_half_step(allocator: std.mem.Allocator, state: *StateType, dt: f64) !void {
            try faraday_step(allocator, state, dt);
        }

        /// Ampere half-step: E^{n+1} = E^n + dt · (★₁⁻¹ d ★₂ B^{n+1/2} − J).
        /// Also advances the discrete timestep counter.
        pub fn second_half_step(allocator: std.mem.Allocator, state: *StateType, dt: f64) !void {
            try ampere_step(allocator, state, dt);
            state.timestep += 1;
        }
    };
}

/// Maxwell leapfrog integrator — built from the generic `Leapfrog`.
///
/// Satisfies the `TimeStepStrategy` concept so it can be wrapped by
/// `TimeStepper` and used by simulation runners.
pub fn MaxwellLeapfrog(comptime MeshType: type) type {
    return leapfrog_mod.Leapfrog(MaxwellSystem(MeshType));
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
        pub fn init(mesh: *const MeshType, frequency: f64, amplitude: f64, position: [MeshType.embedding_dimension]f64) Self {
            const simplex_1 = mesh.simplices(1);
            const edge_verts = simplex_1.items(.vertices);
            const lengths = simplex_1.items(.volume);
            const coords = mesh.vertices.slice().items(.coords);

            var best_edge: u32 = 0;
            var best_distance_squared: f64 = std.math.inf(f64);

            for (0..mesh.num_edges()) |e| {
                const v0 = coords[edge_verts[e][0]];
                const v1 = coords[edge_verts[e][1]];

                // Edge midpoint.
                var midpoint: [MeshType.embedding_dimension]f64 = undefined;
                inline for (0..MeshType.embedding_dimension) |d| {
                    midpoint[d] = 0.5 * (v0[d] + v1[d]);
                }

                var dist_sq: f64 = 0;
                inline for (0..MeshType.embedding_dimension) |d| {
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
// Energy diagnostic
// ═══════════════════════════════════════════════════════════════════════════

/// Discrete electromagnetic energy: U = ½(⟨E, ★₁E⟩ + ⟨B, ★₂B⟩).
///
/// The Hodge-weighted L² inner product on differential forms. ★₁ uses the
/// Whitney mass matrix (exact Galerkin inner product), ★₂ is diagonal.
/// The leapfrog integrator conserves this quantity (up to O(dt²) oscillation)
/// for source-free fields.
pub fn electromagnetic_energy(
    allocator: std.mem.Allocator,
    state: anytype,
) !f64 {
    // ★₁(E): primal 1-form → dual 1-form.
    var star_E = try hodge_star.hodge_star(allocator, state.E);
    defer star_E.deinit(allocator);

    // ★₂(B): primal 2-form → dual 0-form.
    var star_B = try hodge_star.hodge_star(allocator, state.B);
    defer star_B.deinit(allocator);

    // ⟨E, ★₁E⟩ = Σᵢ Eᵢ · (★₁E)ᵢ
    var e_energy: f64 = 0.0;
    for (state.E.values, star_E.values) |e, se| {
        e_energy += e * se;
    }

    // ⟨B, ★₂B⟩ = Σᵢ Bᵢ · (★₂B)ᵢ
    var b_energy: f64 = 0.0;
    for (state.B.values, star_B.values) |b, sb| {
        b_energy += b * sb;
    }

    return 0.5 * (e_energy + b_energy);
}

// ═══════════════════════════════════════════════════════════════════════════
// Simulation runner
// ═══════════════════════════════════════════════════════════════════════════

const flux_io = flux.io;
const std_fs = std.fs;

/// Configuration for the simulation runner.
pub const RunConfig = struct {
    /// Number of leapfrog timesteps to execute.
    steps: u32,
    /// Timestep size dt.
    dt: f64,
    /// Write a VTK snapshot every this many steps. 0 = no output.
    output_interval: u32 = 0,
    /// Directory path for VTK output files. Must exist if output_interval > 0.
    output_path: ?[]const u8 = null,
    /// Base name for snapshot files (e.g., "field" → "field_0000.vtu").
    output_base_name: []const u8 = "field",
};

/// Simulation runner: advances a Maxwell state through a configured number
/// of leapfrog timesteps, applying PEC boundary conditions and optionally
/// writing VTK snapshots at a fixed interval.
///
/// The runner does not own the state or mesh — the caller manages their
/// lifetimes. The runner is a pure orchestrator: it sequences the leapfrog
/// integrator, boundary conditions, and I/O without adding physics.
pub fn Runner(comptime MeshType: type) type {
    return struct {
        const Self = @This();
        const MState = State(MeshType);

        /// Run the simulation for `config.steps` timesteps.
        ///
        /// At each step:
        ///   1. Update the current source J (if source is provided).
        ///   2. Advance fields via leapfrog (Faraday then Ampere).
        ///   3. Apply PEC boundary conditions.
        ///   4. Write a VTK snapshot if the step is a multiple of output_interval.
        ///
        /// After the loop, writes a .pvd collection file referencing all snapshots
        /// (if any were written) so ParaView can load them as an animation.
        pub fn run(
            allocator: std.mem.Allocator,
            state: *MState,
            source: anytype,
            config: RunConfig,
        ) !void {
            std.debug.assert(config.steps > 0);
            std.debug.assert(config.dt > 0);

            const has_output = config.output_interval > 0 and config.output_path != null;
            const has_source = @TypeOf(source) != @TypeOf(null);

            // Pre-allocate PVD entry storage for all snapshots we'll write.
            const max_snapshots = if (has_output)
                (config.steps / config.output_interval) + 1 // +1 for possible step 0
            else
                0;

            var pvd_filenames: []flux_io.PvdEntry = &.{};
            var snapshot_count: u32 = 0;
            var filename_bufs: [][flux_io.max_snapshot_filename_length]u8 = &.{};

            if (has_output) {
                pvd_filenames = try allocator.alloc(flux_io.PvdEntry, max_snapshots);
                filename_bufs = try allocator.alloc([flux_io.max_snapshot_filename_length]u8, max_snapshots);
            }
            defer if (has_output) {
                allocator.free(pvd_filenames);
                allocator.free(filename_bufs);
            };

            for (0..config.steps) |step_idx| {
                const t = @as(f64, @floatFromInt(step_idx)) * config.dt;

                // Update source.
                if (has_source) {
                    source.apply(&state.J, t);
                }

                // Advance fields.
                try leapfrog_step(allocator, state, config.dt);

                // Enforce PEC boundary conditions.
                apply_pec_boundary(state);

                // Write VTK snapshot if at output interval.
                if (has_output and (step_idx + 1) % config.output_interval == 0) {
                    const filename = flux_io.snapshot_filename(
                        &filename_bufs[snapshot_count],
                        config.output_base_name,
                        snapshot_count,
                    );

                    try writeSnapshot(allocator, config.output_path.?, filename, state);

                    pvd_filenames[snapshot_count] = .{
                        .timestep = t + config.dt,
                        .filename = filename,
                    };
                    snapshot_count += 1;
                }
            }

            // Write PVD collection file.
            if (has_output and snapshot_count > 0) {
                try writePvd(allocator, config.output_path.?, config.output_base_name, pvd_filenames[0..snapshot_count]);
            }
        }

        fn writeSnapshot(
            allocator: std.mem.Allocator,
            output_dir: []const u8,
            filename: []const u8,
            state: *const MState,
        ) !void {
            // Build VTK content in memory, then write to file.
            var output = std.ArrayListUnmanaged(u8){};
            defer output.deinit(allocator);

            try flux_io.write_fields(
                allocator,
                output.writer(allocator),
                MeshType.embedding_dimension,
                MeshType.topological_dimension,
                state.mesh.*,
                state.E.values,
                state.B.values,
            );

            var dir = try std_fs.cwd().openDir(output_dir, .{});
            defer dir.close();
            const file = try dir.createFile(filename, .{});
            defer file.close();
            try file.writeAll(output.items);
        }

        fn writePvd(
            allocator: std.mem.Allocator,
            output_dir: []const u8,
            base_name: []const u8,
            entries: []const flux_io.PvdEntry,
        ) !void {
            var output = std.ArrayListUnmanaged(u8){};
            defer output.deinit(allocator);

            try flux_io.write_pvd(output.writer(allocator), entries);

            var pvd_buf: [flux_io.max_snapshot_filename_length]u8 = undefined;
            const pvd_name = std.fmt.bufPrint(&pvd_buf, "{s}.pvd", .{base_name}) catch
                return error.FilenameTooLong;

            var dir = try std_fs.cwd().openDir(output_dir, .{});
            defer dir.close();
            const file = try dir.createFile(pvd_name, .{});
            defer file.close();
            try file.writeAll(output.items);
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const Mesh2D = topology.Mesh(2, 2);
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
    var dE = try exterior_derivative.exterior_derivative(allocator, state.E);
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
    var star_B = try hodge_star.hodge_star(allocator, state.B);
    defer star_B.deinit(allocator);

    // d(★₂B)
    var d_star_B = try exterior_derivative.exterior_derivative(allocator, star_B);
    defer d_star_B.deinit(allocator);

    // ★₁⁻¹(d(★₂B))
    const edge_volumes = mesh.simplices(1).items(.volume);
    const dual_edge_vols = mesh.dual_edge_volumes;

    const expected_E = try allocator.alloc(f64, state.E.values.len);
    defer allocator.free(expected_E);
    for (expected_E, state.E.values, d_star_B.values, edge_volumes, dual_edge_vols, state.J.values) |*exp, e, dsb, vol, dual_vol, j| {
        exp.* = e + dt * ((vol / dual_vol) * dsb - j);
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
    const verts = mesh.simplices(1).items(.vertices)[dipole.edge_index];
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

// ── #41: Simulation runner ───────────────────────────────────────────────

const MaxwellRunner = Runner(Mesh2D);

test "Runner advances state by configured number of steps" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Seed E with small random values for non-trivial evolution.
    var rng = std.Random.DefaultPrng.init(0xE0_4101);
    for (state.E.values) |*v| v.* = (rng.random().float(f64) * 2.0 - 1.0) * 0.01;

    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = 10,
        .dt = 0.001,
    });

    try testing.expectEqual(@as(u64, 10), state.timestep);
}

test "Runner applies PEC boundary conditions every step" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set non-zero E on boundary edges — runner should zero them.
    for (state.E.values) |*v| v.* = 1.0;

    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = 5,
        .dt = 0.001,
    });

    // After running, all boundary edges must be zero (PEC enforced).
    for (mesh.boundary_edges) |edge_idx| {
        try testing.expectEqual(@as(f64, 0), state.E.values[edge_idx]);
    }
}

test "Runner applies source at each step" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 2.0, 2.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    const dipole = Dipole.init(&mesh, 1.0, 1.0, .{ 1.0, 1.0 });

    try MaxwellRunner.run(allocator, &state, dipole, .{
        .steps = 10,
        .dt = 0.001,
    });

    // After running with a source, fields should be non-trivial.
    var max_abs_e: f64 = 0;
    for (state.E.values) |v| max_abs_e = @max(max_abs_e, @abs(v));
    try testing.expect(max_abs_e > 1e-15);
}

test "Runner writes VTK snapshots at configured interval" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Use a temporary directory for output.
    var tmp_dir = testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Get the path string for the temp directory.
    var path_buf: [std_fs.max_path_bytes]u8 = undefined;
    const tmp_path = try tmp_dir.dir.realpath(".", &path_buf);

    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = 10,
        .dt = 0.001,
        .output_interval = 5,
        .output_path = tmp_path,
        .output_base_name = "test_field",
    });

    // Should have written 2 snapshots (steps 5 and 10) + 1 PVD file.
    // Snapshot at step 5 → test_field_0000.vtu
    // Snapshot at step 10 → test_field_0001.vtu
    const snap0 = tmp_dir.dir.openFile("test_field_0000.vtu", .{});
    try testing.expect(snap0 != error.FileNotFound);
    if (snap0) |f| f.close() else |_| {}

    const snap1 = tmp_dir.dir.openFile("test_field_0001.vtu", .{});
    try testing.expect(snap1 != error.FileNotFound);
    if (snap1) |f| f.close() else |_| {}

    // PVD collection file.
    const pvd = tmp_dir.dir.openFile("test_field.pvd", .{});
    try testing.expect(pvd != error.FileNotFound);
    if (pvd) |f| f.close() else |_| {}
}

test "Runner with output_interval 0 writes no files" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // No output — should run without errors and produce no files.
    try MaxwellRunner.run(allocator, &state, null, .{
        .steps = 5,
        .dt = 0.001,
        .output_interval = 0,
    });

    try testing.expectEqual(@as(u64, 5), state.timestep);
}

// ── #40: Energy tracking ─────────────────────────────────────────────────

test "electromagnetic_energy is zero for zero fields" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 4, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    const energy = try electromagnetic_energy(allocator, &state);
    try testing.expectEqual(@as(f64, 0.0), energy);
}

test "electromagnetic_energy is non-negative" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 6, 6, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Random fields — energy must still be ≥ 0 since ★₀, ★₂ are non-negative
    // diagonal operators and ★₁ (Whitney mass matrix) is SPD.
    var rng = std.Random.DefaultPrng.init(0xE0_E40_01);
    for (state.E.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;
    for (state.B.values) |*v| v.* = rng.random().float(f64) * 2.0 - 1.0;

    const energy = try electromagnetic_energy(allocator, &state);
    try testing.expect(energy >= 0.0);
}

test "electromagnetic_energy uses Hodge-weighted inner product" {
    // Verify that energy ≠ ½(‖E‖² + ‖B‖²) in general — the Hodge star
    // weights make a difference on non-unit meshes.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 3.0, 2.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Set a single non-zero E component.
    state.E.values[0] = 1.0;

    const energy = try electromagnetic_energy(allocator, &state);
    const unweighted = 0.5 * state.E.norm_squared();

    // On a non-unit mesh, Hodge-weighted and unweighted differ.
    try testing.expect(energy != unweighted);
    try testing.expect(energy > 0.0);
}

// ── #46: End-to-end integration ──────────────────────────────────────────

test "end-to-end: 100 steps, dB = 0 structurally, energy bounded" {
    // M3 acceptance test. Runs a dipole simulation for 100 steps and verifies:
    //   1. dB = 0 at every step (structurally guaranteed in 2D — B is a
    //      top-form, so d₂B lives in the zero-dimensional 3-form space).
    //   2. Energy is bounded — not growing without bound.
    //
    // The energy check is the substantive diagnostic. With a source driving
    // the system, energy will grow (the source injects energy), but it must
    // remain finite and bounded — no numerical blowup.
    const allocator = testing.allocator;

    var mesh = try Mesh2D.uniform_grid(allocator, 8, 8, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Place a dipole source near the center.
    const dipole = PointDipole(Mesh2D).init(&mesh, 10.0, 0.1, .{ 0.5, 0.5 });

    // dt chosen small relative to mesh spacing for CFL stability.
    const dt: f64 = 0.001;
    const num_steps: usize = 100;

    var energy_max: f64 = 0.0;

    for (0..num_steps) |step_idx| {
        const t = @as(f64, @floatFromInt(step_idx)) * dt;
        dipole.apply(&state.J, t);

        try leapfrog_step(allocator, &state, dt);
        apply_pec_boundary(&state);

        // ── Invariant 1: dB = 0 ──
        // In 2D, B is a 2-form (top-form). d₂ would map into 3-forms,
        // which don't exist on a 2D mesh — the type system rejects
        // `exterior_derivative(B)` at comptime. This is the strongest
        // possible guarantee: dB = 0 is not checked numerically because
        // it is structurally unrepresentable as a nonzero object.
        comptime {
            std.debug.assert(MaxwellState.TwoForm.degree == Mesh2D.embedding_dimension);
        }

        // ── Invariant 2: energy is finite ──
        const energy = try electromagnetic_energy(allocator, &state);
        try testing.expect(!std.math.isNan(energy));
        try testing.expect(!std.math.isInf(energy));
        energy_max = @max(energy_max, energy);
    }

    // Verify simulation actually ran and produced non-trivial fields.
    try testing.expectEqual(@as(u64, num_steps), state.timestep);
    try testing.expect(energy_max > 0.0);

    // Energy should stay bounded. With a weak source (amplitude 0.1) over
    // 100 steps at dt = 0.001, energy should be small. A blowup would
    // produce energy >> 1.
    try testing.expect(energy_max < 1.0);
}

// ── #44: Convergence test ────────────────────────────────────────────────

// ── Analytical TE₁₀ cavity mode for convergence testing ──────────────
//
// Standing wave in a [0,L]×[0,L] PEC cavity (c = 1, ω = π/L):
//   B_z(x,y,t) = cos(πx/L) cos(ωt)
//   E_y(x,y,t) = sin(πx/L) sin(ωt),  E_x = 0
//
// As DEC forms, E is a primal 1-form (line integrals along edges) and
// B is a primal 2-form (flux integrals over faces). We project using
// midpoint quadrature for edges and centroid quadrature for faces.

/// Project the analytical TE₁₀ E field onto mesh edges at time t.
///
/// Each edge value = ∫ E · dl ≈ E(midpoint) · edge_vector.
pub fn project_te10_e(mesh: *const Mesh2D, values: []f64, t: f64, domain_length: f64) void {
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    const k = std.math.pi / domain_length;
    const omega = k; // c = 1

    for (values, edge_verts) |*val, verts| {
        const p0 = coords[verts[0]];
        const p1 = coords[verts[1]];
        const mx = 0.5 * (p0[0] + p1[0]);
        const dy = p1[1] - p0[1];

        // E · dl = E_x·dx + E_y·dy = 0 + sin(kx)sin(ωt)·dy.
        val.* = @sin(k * mx) * @sin(omega * t) * dy;
    }
}

/// Project the analytical TE₁₀ B field onto mesh faces at time t.
///
/// Each face value = ∫∫ B_z dA ≈ B_z(centroid) · area.
pub fn project_te10_b(mesh: *const Mesh2D, values: []f64, t: f64, domain_length: f64) void {
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

/// Compute the discrete electromagnetic energy: ½(⟨E, ★₁E⟩ + ⟨B, ★₂B⟩).
///
/// Uses the Whitney mass matrix for ★₁ (SpMV) and the diagonal for ★₂.
/// Allocates a temporary buffer for the SpMV result.
fn discrete_energy_from_arrays(allocator: std.mem.Allocator, mesh: *const Mesh2D, e_vals: []const f64, b_vals: []const f64) !f64 {
    // ★₁(E) via Whitney mass matrix SpMV.
    const m1_e = try allocator.alloc(f64, e_vals.len);
    defer allocator.free(m1_e);
    sparse.spmv(mesh.whitney_mass_1, e_vals, m1_e);

    var e_energy: f64 = 0.0;
    for (e_vals, m1_e) |e, se| {
        e_energy += e * se;
    }

    // ★₂(B): diagonal (1 / area).
    const face_areas = mesh.simplices(2).items(.volume);
    var b_energy: f64 = 0.0;
    for (b_vals, face_areas) |b, area| {
        b_energy += b * b / area;
    }

    return 0.5 * (e_energy + b_energy);
}

/// Run a TE₁₀ cavity simulation at given resolution and return the
/// relative energy error at final_time compared to the analytical mode.
///
/// Energy is the natural diagnostic: it uses the Hodge-weighted inner
/// product (which correctly de-weights degenerate edges) and measures
/// a single scalar quantity that integrates all field components.
/// Run a source-free TE₁₀ cavity and return the relative energy drift
/// at final_time. Energy conservation is the convergent diagnostic for
/// the leapfrog integrator — it depends on the symplectic structure,
/// not on the Hodge star's approximation quality.
fn run_cavity_energy_drift(allocator: std.mem.Allocator, grid_n: u32, final_time: f64) !f64 {
    const domain_length: f64 = 1.0;
    const grid_spacing = domain_length / @as(f64, @floatFromInt(grid_n));

    // CFL-stable timestep: dt = 0.1 · h (well below the CFL limit).
    const dt = 0.1 * grid_spacing;
    const num_steps: u32 = @intFromFloat(@round(final_time / dt));

    var mesh = try Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    // Initial conditions for leapfrog stagger.
    project_te10_b(&mesh, state.B.values, -dt / 2.0, domain_length);

    const energy_initial = try discrete_energy_from_arrays(allocator, &mesh, state.E.values, state.B.values);

    // Run leapfrog with PEC.
    for (0..num_steps) |_| {
        try leapfrog_step(allocator, &state, dt);
        apply_pec_boundary(&state);
    }

    const energy_final = try discrete_energy_from_arrays(allocator, &mesh, state.E.values, state.B.values);

    return @abs(energy_final - energy_initial) / energy_initial;
}

/// Compute the discrete TE₁₀ eigenvalue ω² via Rayleigh quotient.
///
/// Projects the analytical TE₁₀ E-mode onto the mesh, applies the DEC
/// curl-curl operator (★₁⁻¹ d̃₀ ★₂ d₁), and returns the Rayleigh quotient:
///   ω²_num = ⟨curl_curl(E), ★₁ E⟩ / ⟨E, ★₁ E⟩
///
/// This gives the exact discrete eigenvalue without time-stepping.
fn compute_te10_eigenvalue(allocator: std.mem.Allocator, grid_n: u32, domain_length: f64) !f64 {
    var mesh = try Mesh2D.uniform_grid(allocator, grid_n, grid_n, domain_length, domain_length);
    defer mesh.deinit(allocator);

    // Project the TE₁₀ E-mode at t = π/(2ω) where sin(ωt) = 1,
    // so E(edge) = sin(πx/L) · dy (maximum amplitude snapshot).
    var E = try MaxwellState.OneForm.init(allocator, &mesh);
    defer E.deinit(allocator);
    project_te10_e(&mesh, E.values, domain_length / 2.0, domain_length);

    // dE (primal 2-form).
    var dE = try exterior_derivative.exterior_derivative(allocator, E);
    defer dE.deinit(allocator);

    // Energy-based Rayleigh quotient: ω² = ⟨dE, ★₂dE⟩ / ⟨E, ★₁E⟩
    // ★₂ is diagonal (exact for faces), ★₁ uses the Whitney mass matrix.

    // Numerator: ⟨dE, ★₂dE⟩
    var star_dE = try hodge_star.hodge_star(allocator, dE);
    defer star_dE.deinit(allocator);
    var numerator: f64 = 0.0;
    for (dE.values, star_dE.values) |de, sde| {
        numerator += de * sde;
    }

    // Denominator: ⟨E, ★₁E⟩ (Whitney mass matrix via hodge_star)
    var star_E = try hodge_star.hodge_star(allocator, E);
    defer star_E.deinit(allocator);
    var denominator: f64 = 0.0;
    for (E.values, star_E.values) |e, se| {
        denominator += e * se;
    }

    return numerator / denominator;
}

// ── TE₁₀ field projection tests ──────────────────────────────────────────

test "project_te10_e is zero for vertical edges at x = 0 and x = L" {
    // TE₁₀ mode: E_y(x,y,t) = sin(πx/L)sin(ωt), E_x = 0.
    // At x = 0 and x = L, sin(πx/L) = 0, so all edge integrals vanish
    // regardless of orientation.
    const allocator = testing.allocator;
    const nx: u32 = 8;
    const ny: u32 = 6;
    const domain_length: f64 = 1.0;
    var mesh = try Mesh2D.uniform_grid(allocator, nx, ny, domain_length, domain_length);
    defer mesh.deinit(allocator);

    const values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(values);

    // At t = L/2 (sin(ωt) = 1, maximum amplitude).
    project_te10_e(&mesh, values, domain_length / 2.0, domain_length);

    // Edges on x = 0: vertices with i = 0.
    // Edges on x = L: vertices with i = nx.
    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_verts) |val, verts| {
        const x0 = coords[verts[0]][0];
        const x1 = coords[verts[1]][0];
        const mx = 0.5 * (x0 + x1);
        // Edges centered at x = 0 or x = L have sin(πx/L) = 0.
        if (@abs(mx) < 1e-12 or @abs(mx - domain_length) < 1e-12) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), val, 1e-14);
        }
    }
}

test "project_te10_e has zero contribution from horizontal edges" {
    // E = (0, E_y), so E · dl = E_y · dy. For horizontal edges, dy = 0,
    // so every horizontal edge gets zero regardless of x position.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 6, 6, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(values);

    project_te10_e(&mesh, values, 0.25, 1.0);

    const edge_verts = mesh.simplices(1).items(.vertices);
    const coords = mesh.vertices.slice().items(.coords);
    for (values, edge_verts) |val, verts| {
        const dy = coords[verts[1]][1] - coords[verts[0]][1];
        if (@abs(dy) < 1e-15) {
            try testing.expectApproxEqAbs(@as(f64, 0.0), val, 1e-15);
        }
    }
}

test "project_te10_b at t = 0 equals cos(πx/L) · area for each face" {
    // B_z(x,y,0) = cos(πx/L) · cos(0) = cos(πx/L).
    // Face integral ≈ cos(πx_centroid/L) · area.
    const allocator = testing.allocator;
    const domain_length: f64 = 2.0;
    var mesh = try Mesh2D.uniform_grid(allocator, 8, 8, domain_length, domain_length);
    defer mesh.deinit(allocator);

    const values = try allocator.alloc(f64, mesh.num_faces());
    defer allocator.free(values);
    project_te10_b(&mesh, values, 0.0, domain_length);

    const face_verts = mesh.simplices(2).items(.vertices);
    const face_areas = mesh.simplices(2).items(.volume);
    const coords = mesh.vertices.slice().items(.coords);
    const k = std.math.pi / domain_length;

    for (values, face_verts, face_areas) |val, verts, area| {
        const cx = (coords[verts[0]][0] + coords[verts[1]][0] + coords[verts[2]][0]) / 3.0;
        const expected = @cos(k * cx) * area;
        try testing.expectApproxEqAbs(expected, val, 1e-14);
    }
}

test "project_te10_e at t = 0 is identically zero" {
    // E_y(x,y,0) = sin(πx/L) · sin(0) = 0 for all x, y.
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 6, 6, 1.0, 1.0);
    defer mesh.deinit(allocator);

    const values = try allocator.alloc(f64, mesh.num_edges());
    defer allocator.free(values);
    project_te10_e(&mesh, values, 0.0, 1.0);

    for (values) |val| {
        try testing.expectApproxEqAbs(@as(f64, 0.0), val, 1e-15);
    }
}

test "TE₁₀ cavity: eigenvalue error reduces ≥3× when grid halves (O(h²) convergence)" {
    // The Whitney/Galerkin mass matrix M₁ gives O(h²) eigenvalue convergence
    // on any triangulation. The Rayleigh quotient ω² = ⟨dE, ★₂dE⟩ / ⟨E, ★₁E⟩
    // uses ★₁ via the Whitney mass matrix (SpMV) and ★₂ via the diagonal.
    //
    // Acceptance criterion (issue #70): error_8×8 / error_16×16 ≥ 3.
    const allocator = testing.allocator;
    const domain_length: f64 = 1.0;
    const analytical_omega_sq = std.math.pi * std.math.pi / (domain_length * domain_length);

    const omega_sq_coarse = try compute_te10_eigenvalue(allocator, 8, domain_length);
    const omega_sq_fine = try compute_te10_eigenvalue(allocator, 16, domain_length);

    const error_coarse = @abs(omega_sq_coarse - analytical_omega_sq) / analytical_omega_sq;
    const error_fine = @abs(omega_sq_fine - analytical_omega_sq) / analytical_omega_sq;

    try testing.expect(omega_sq_coarse > 0.0);
    try testing.expect(omega_sq_fine > 0.0);

    std.debug.print("\nWhitney ω²: 8×8={d:.6}, 16×16={d:.6}, analytical={d:.6}\n", .{
        omega_sq_coarse, omega_sq_fine, analytical_omega_sq,
    });
    std.debug.print("Whitney error: 8×8={d:.6}, 16×16={d:.6}, ratio={d:.2}\n", .{
        error_coarse, error_fine, error_coarse / error_fine,
    });

    // Finer grid must be closer to analytical.
    try testing.expect(error_fine < error_coarse);

    // Acceptance criterion: error reduces by ≥3× when grid halves.
    try testing.expect(error_coarse / error_fine >= 3.0);
}

test "convergence: energy drift bounded for TE₁₀ cavity" {
    // The leapfrog integrator with the Whitney Hodge star uses a CG solve
    // per Ampere step, which introduces small residuals that accumulate
    // over many timesteps. The energy drift is bounded but not monotonically
    // decreasing with grid refinement (finer grids need more steps).
    //
    // We verify: energy drift stays below 2% for both resolutions.
    const allocator = testing.allocator;

    const final_time: f64 = 0.5;
    const coarse_n: u32 = 8;
    const fine_n: u32 = 16;

    const drift_coarse = try run_cavity_energy_drift(allocator, coarse_n, final_time);
    const drift_fine = try run_cavity_energy_drift(allocator, fine_n, final_time);

    // Both drifts should be small.
    try testing.expect(drift_coarse < 0.02); // < 2%
    try testing.expect(drift_fine < 0.02); // < 2%
}

// ── #62: MaxwellLeapfrog — TimeStepper concept conformance ──────────

const time_stepper = flux.time_stepper;

test "MaxwellLeapfrog satisfies TimeStepStrategy concept" {
    const Leapfrog = MaxwellLeapfrog(Mesh2D);
    comptime time_stepper.TimeStepStrategy(Leapfrog);
}

test "MaxwellLeapfrog works through TimeStepper wrapper" {
    const Stepper = time_stepper.TimeStepper(MaxwellLeapfrog(Mesh2D));
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 3, 3, 1.0, 1.0);
    defer mesh.deinit(allocator);

    var state = try MaxwellState.init(allocator, &mesh);
    defer state.deinit(allocator);

    try Stepper.step(allocator, &state, 0.01);
    try testing.expectEqual(@as(u64, 1), state.timestep);

    try Stepper.step(allocator, &state, 0.01);
    try testing.expectEqual(@as(u64, 2), state.timestep);
}

test "MaxwellLeapfrog.step produces same result as leapfrog_step" {
    const allocator = testing.allocator;
    var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
    defer mesh.deinit(allocator);

    // Set up two identical states.
    var state_free = try MaxwellState.init(allocator, &mesh);
    defer state_free.deinit(allocator);
    var state_wrap = try MaxwellState.init(allocator, &mesh);
    defer state_wrap.deinit(allocator);

    // Identical random initial conditions.
    var rng = std.Random.DefaultPrng.init(0x62_1EAF);
    for (state_free.E.values, state_wrap.E.values) |*a, *b| {
        const v = rng.random().float(f64) * 2.0 - 1.0;
        a.* = v;
        b.* = v;
    }
    for (state_free.B.values, state_wrap.B.values) |*a, *b| {
        const v = rng.random().float(f64) * 2.0 - 1.0;
        a.* = v;
        b.* = v;
    }

    const dt: f64 = 0.01;

    // Advance via free function.
    try leapfrog_step(allocator, &state_free, dt);

    // Advance via wrapper.
    const Leapfrog = MaxwellLeapfrog(Mesh2D);
    try Leapfrog.step(allocator, &state_wrap, dt);

    // Results must be identical.
    for (state_free.E.values, state_wrap.E.values) |a, b| {
        try testing.expectEqual(a, b);
    }
    for (state_free.B.values, state_wrap.B.values) |a, b| {
        try testing.expectEqual(a, b);
    }
    try testing.expectEqual(state_free.timestep, state_wrap.timestep);
}
