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
//! var mesh = try Mesh2D.plane(allocator, 4, 3, 2.0, 1.5);
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

pub const forms = @import("forms/root.zig");
pub const io = @import("io/root.zig");
pub const math = struct {
    pub const sparse = @import("math/sparse.zig");
    pub const cg = @import("math/cg.zig");
    pub const linear_system = @import("math/linear_system.zig");
};
const laplacian_mod = @import("operators/laplacian.zig");
pub const operators = struct {
    pub const context = @import("operators/context.zig");
    pub const dec = struct {
        pub const exterior_derivative = @import("operators/exterior_derivative.zig");
        pub const laplacian = laplacian_mod.dec;
    };
    pub const feec = struct {
        pub const codifferential = @import("operators/codifferential.zig");
        pub const hodge_star = @import("operators/hodge_star.zig");
        pub const laplacian = laplacian_mod.feec;
        pub const weak_form = @import("operators/weak_form.zig");
        pub const whitney_mass = @import("operators/whitney_mass.zig");
    };
    pub const bridges = @import("operators/bridges.zig");
    pub const boundary_conditions = @import("operators/boundary_conditions.zig");
    pub const codifferential = @import("operators/codifferential.zig");
    pub const compose = @import("operators/compose.zig");
    pub const exterior_derivative = @import("operators/exterior_derivative.zig");
    pub const hodge_star = @import("operators/hodge_star.zig");
    pub const observers = @import("operators/observers.zig");
    pub const poisson = @import("operators/poisson.zig");
    pub const weak_form = @import("operators/weak_form.zig");
    pub const whitney_mass = @import("operators/whitney_mass.zig");
    pub const wedge_product = @import("operators/wedge_product.zig");
};
pub const topology = @import("topology/mesh.zig");
pub const evolution = @import("evolution/root.zig");
pub const listeners = evolution.listeners;
pub const integrators = struct {
    pub const Leapfrog = @import("integrators/leapfrog.zig").Leapfrog;
    pub const ForwardEuler = @import("integrators/forward_euler.zig").ForwardEuler;
    pub const leapfrog = @import("integrators/leapfrog.zig");
    pub const forward_euler = @import("integrators/forward_euler.zig");
};
pub const concepts = struct {
    pub const mesh = @import("concepts/mesh.zig");
};

test "public API exposes explicit DEC and FEEC operator families" {
    const testing = @import("std").testing;

    try testing.expect(@hasDecl(@This().operators, "dec"));
    try testing.expect(@hasDecl(@This().operators, "feec"));
    try testing.expect(@hasDecl(@This().operators, "context"));
    try testing.expect(@hasDecl(@This().operators.dec, "laplacian"));
    try testing.expect(@hasDecl(@This().operators.feec, "laplacian"));
    try testing.expect(!@hasDecl(@This().operators, "laplacian"));
}

test "public API exposes evolution as the execution root module" {
    const testing = @import("std").testing;

    try testing.expect(@hasDecl(@This(), "evolution"));
    try testing.expect(@hasDecl(@This(), "listeners"));
    try testing.expect(!@hasDecl(@This().integrators, "evolution"));
}

test "operator context owns assembled operators directly" {
    const testing = @import("std").testing;
    const Mesh2D = topology.Mesh(2, 2);
    const OperatorContext = @This().operators.context.OperatorContext(Mesh2D);

    try testing.expect(@hasDecl(OperatorContext, "exteriorDerivative"));
    try testing.expect(@hasDecl(OperatorContext, "hodgeStar"));
    try testing.expect(@hasDecl(OperatorContext, "laplacian"));
}

test "forms API exposes FEEC spaces while keeping Cochain storage-only" {
    const testing = @import("std").testing;
    const Mesh2D = topology.Mesh(2, 2);
    const C1 = @This().forms.Cochain(Mesh2D, 1, @This().forms.Primal);

    try testing.expect(@hasDecl(@This().forms, "feec"));
    try testing.expect(@hasDecl(@This().forms.feec, "WhitneySpace"));
    try testing.expect(!@hasDecl(C1, "interpolate"));
    try testing.expect(!@hasDecl(C1, "project"));
    try testing.expect(!@hasDecl(C1, "space"));
}

test "operators API exposes explicit DEC/FEEC bridge operators" {
    const testing = @import("std").testing;

    try testing.expect(@hasDecl(@This().operators.bridges, "WhitneyInterpolation"));
    try testing.expect(@hasDecl(@This().operators.bridges, "DeRhamProjection"));
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
