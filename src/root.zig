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
//! const Mesh2D = flux.topology.Mesh(2, 2);
//! var mesh = try Mesh2D.uniform_grid(allocator, 4, 3, 2.0, 1.5);
//! const operators = try flux.operators.context.OperatorContext(Mesh2D).init(allocator, &mesh);
//! defer operators.deinit();
//!
//! // Create a 0-cochain (scalar field on vertices)
//! var omega = try flux.forms.Cochain(Mesh2D, 0, flux.forms.Primal).init(allocator, &mesh);
//!
//! // Apply operators with compile-time type safety
//! var d_omega = try (try operators.exteriorDerivative(flux.forms.Primal, 0)).apply(allocator, omega);
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
    /// A tetrahedron has zero volume — 3D dual geometry and boundary
    /// operators cannot be computed consistently.
    DegenerateTetrahedron,
    /// An edge is adjacent to zero faces — the mesh is disconnected or
    /// malformed.
    NonManifoldEdge,
    /// A face is adjacent to more than two tetrahedra, so the simplicial
    /// complex is not a manifold-with-boundary.
    NonManifoldFace,
};

// ── Submodule re-exports (for namespaced access) ────────────────────────

pub const forms = @import("forms/cochain.zig");
pub const io = @import("io/root.zig");
pub const math = struct {
    pub const sparse = @import("math/sparse.zig");
    pub const cg = @import("math/cg.zig");
};
pub const operators = struct {
    pub const boundary_conditions = @import("operators/boundary_conditions.zig");
    pub const codifferential = @import("operators/codifferential.zig");
    pub const compose = @import("operators/compose.zig");
    pub const context = @import("operators/context.zig");
    pub const exterior_derivative = @import("operators/exterior_derivative.zig");
    pub const hodge_star = @import("operators/hodge_star.zig");
    pub const implicit_system = @import("operators/implicit_system.zig");
    pub const laplacian = @import("operators/laplacian.zig");
    pub const observers = @import("operators/observers.zig");
    pub const poisson = @import("operators/poisson.zig");
    pub const whitney_mass = @import("operators/whitney_mass.zig");
    pub const wedge_product = @import("operators/wedge_product.zig");
};
pub const topology = @import("topology/mesh.zig");
pub const time_stepper = @import("time_stepper.zig");
pub const integrators = struct {
    pub const evolution = @import("integrators/evolution.zig");
    pub const leapfrog = @import("integrators/leapfrog.zig");
    pub const forward_euler = @import("integrators/forward_euler.zig");
};
pub const concepts = struct {
    pub const mesh = @import("concepts/mesh.zig");
};

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
