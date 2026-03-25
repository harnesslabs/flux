//! flux — composable, type-safe PDE solver framework.
//!
//! flux assembles high-fidelity numerical solvers by composing discrete
//! operators on function spaces defined over simplicial meshes. The central
//! design principle is that incompatible operator compositions (e.g., applying
//! a 1-form operator to a 2-form) are caught at compile time, not runtime.
//!
//! ## Quick start
//!
//! ```zig
//! const flux = @import("flux");
//!
//! // Build a mesh
//! var mesh = try flux.Mesh(2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
//!
//! // Create a 0-cochain (scalar field on vertices)
//! var omega = try flux.Cochain(flux.Mesh(2), 0, flux.Primal).init(allocator, &mesh);
//!
//! // Apply operators with compile-time type safety
//! var d_omega = try flux.exterior_derivative(allocator, omega);
//! var star_d  = try flux.hodge_star(allocator, d_omega);
//! ```
//!
//! ## Modules
//!
//! - `topology` — simplicial meshes, boundary operators, geometric data
//! - `forms` — discrete k-forms (cochains) parameterized on mesh and degree
//! - `operators` — d, ★, ★⁻¹, Δ, and composition API
//! - `math` — sparse linear algebra primitives (CSR matrices)
//! - `io` — VTK `.vtu` serialization for ParaView visualization

// ── Submodule re-exports (for namespaced access) ────────────────────────

pub const forms = @import("forms/cochain.zig");
pub const io = @import("io/vtk.zig");
pub const math = @import("math/sparse.zig");
pub const operators = struct {
    pub const compose = @import("operators/compose.zig");
    pub const exterior_derivative = @import("operators/exterior_derivative.zig");
    pub const hodge_star = @import("operators/hodge_star.zig");
    pub const laplacian = @import("operators/laplacian.zig");
};
pub const topology = @import("topology/mesh.zig");
pub const em = @import("em/maxwell.zig");

/// Electromagnetic field state (E, B, J) on a simplicial mesh.
pub const MaxwellState = em.State;

/// Apply PEC boundary conditions — zero E on all boundary edges.
pub const apply_pec_boundary = em.apply_pec_boundary;

/// Point dipole current source — sinusoidal J on the nearest edge.
pub const PointDipole = em.PointDipole;

/// Simulation runner — drives the leapfrog loop with configurable output.
pub const Runner = em.Runner;

/// Configuration for the simulation runner.
pub const RunConfig = em.RunConfig;

/// Discrete electromagnetic energy: ½(⟨E, ★₁E⟩ + ⟨B, ★₂B⟩).
pub const electromagnetic_energy = em.electromagnetic_energy;

/// Write E and B fields to a VTK .vtu file (E projected onto faces).
pub const write_fields = io.write_fields;

/// Project a 1-form (edge values) to per-face averages for VTK export.
pub const project_edges_to_faces = io.project_edges_to_faces;

// ── Top-level convenience re-exports ────────────────────────────────────
// These allow `const flux = @import("flux"); flux.Cochain(...)` without
// navigating the module hierarchy.

/// Discrete k-form (cochain) on a simplicial mesh — assigns one real value
/// per k-cell. Parameterized on mesh type, degree k, and duality (Primal/Dual).
pub const Cochain = forms.Cochain;

/// Marker: cochain lives on the primal complex.
pub const Primal = forms.Primal;

/// Marker: cochain lives on the dual complex.
pub const Dual = forms.Dual;

/// Simplicial mesh parameterized on topological dimension.
pub const Mesh = topology.Mesh;

/// Exterior derivative dₖ: Ωᵏ → Ωᵏ⁺¹. Maps k-cochains to (k+1)-cochains
/// via the coboundary operator. Works on both primal and dual cochains.
pub const exterior_derivative = operators.exterior_derivative.exterior_derivative;

/// Hodge star ★ₖ: primal Ωᵏ → dual Ωⁿ⁻ᵏ. Diagonal operator encoding the
/// mesh metric — scales each coefficient by the ratio of dual to primal
/// cell volumes.
pub const hodge_star = operators.hodge_star.hodge_star;

/// Inverse Hodge star ★⁻¹: dual Ωⁿ⁻ᵏ → primal Ωᵏ. Element-wise reciprocal
/// of the Hodge star diagonal.
pub const hodge_star_inverse = operators.hodge_star.hodge_star_inverse;

/// Hodge Laplacian Δₖ = dδ + δd on primal k-cochains. Self-adjoint,
/// positive-semidefinite on 0-forms.
pub const laplacian = operators.laplacian.laplacian;

/// Apply a sequence of DEC operators with automatic intermediate allocation
/// management. Degree/duality mismatches are compile errors.
pub const chain = operators.compose.chain;

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
