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
//! var mesh = try flux.Mesh(2, 2).uniform_grid(allocator, 4, 3, 2.0, 1.5);
//!
//! // Create a 0-cochain (scalar field on vertices)
//! var omega = try flux.Cochain(flux.Mesh(2, 2), 0, flux.Primal).init(allocator, &mesh);
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

// ── Library error set ────────────────────────────────────────────────────

/// Errors arising from mesh data that violates geometric or topological
/// invariants. These indicate problems with input data (e.g., a loaded
/// mesh file), not programmer bugs — programmer bugs are assertions that
/// panic immediately.
pub const Error = error{
    /// A triangle has zero area — cotangent weights and dual geometry
    /// cannot be computed.
    DegenerateTriangle,
    /// An edge is adjacent to zero faces — the mesh is disconnected or
    /// malformed.
    NonManifoldEdge,
};

// ── Submodule re-exports (for namespaced access) ────────────────────────

pub const forms = @import("forms/cochain.zig");
pub const io = @import("io/vtk.zig");
pub const math = struct {
    pub const sparse = @import("math/sparse.zig");
    pub const cg = @import("math/cg.zig");
};
pub const operators = struct {
    pub const compose = @import("operators/compose.zig");
    pub const exterior_derivative = @import("operators/exterior_derivative.zig");
    pub const hodge_star = @import("operators/hodge_star.zig");
    pub const laplacian = @import("operators/laplacian.zig");
    pub const whitney_mass = @import("operators/whitney_mass.zig");
};
pub const topology = @import("topology/mesh.zig");
pub const time_stepper = @import("time_stepper.zig");
pub const integrators = struct {
    pub const leapfrog = @import("integrators/leapfrog.zig");
    pub const forward_euler = @import("integrators/forward_euler.zig");
};
pub const concepts = struct {
    pub const mesh = @import("concepts/mesh.zig");
};

// ── Comptime concepts ───────────────────────────────────────────────────
// Validators that enforce compile-time contracts on user-defined types.
// Each is a function called at comptime that produces @compileError on
// non-conformance.

/// Comptime concept: validate that a type satisfies the TimeStepStrategy contract.
pub const TimeStepStrategy = time_stepper.TimeStepStrategy;

/// Generic time integrator — wraps any conforming TimeStepStrategy.
pub const TimeStepper = time_stepper.TimeStepper;

/// Comptime concept: validate that a type satisfies the Mesh contract.
/// Requires dimension, topological_dimension, entity count accessors, and boundary().
pub const MeshConcept = concepts.mesh.MeshConcept;

/// Generic leapfrog integrator — composes two symplectic half-steps.
pub const Leapfrog = integrators.leapfrog.Leapfrog;

/// Generic forward Euler integrator — single explicit step.
pub const ForwardEuler = integrators.forward_euler.ForwardEuler;

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

/// Simplicial mesh parameterized on `embedding_dimension` and `topological_dimension`.
pub const Mesh = topology.Mesh;

/// Exterior derivative dₖ: Ωᵏ → Ωᵏ⁺¹. Maps k-cochains to (k+1)-cochains
/// via the coboundary operator. Works on both primal and dual cochains.
pub const exterior_derivative = operators.exterior_derivative.exterior_derivative;

/// Hodge star ★ₖ: primal Ωᵏ → dual Ωⁿ⁻ᵏ. Diagonal for k=0,n; Whitney
/// mass matrix (SpMV) for 0 < k < n.
pub const hodge_star = operators.hodge_star.hodge_star;

/// Inverse Hodge star ★⁻¹: dual Ωⁿ⁻ᵏ → primal Ωᵏ. Diagonal for k=0,n;
/// preconditioned CG solve for 0 < k < n.
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
