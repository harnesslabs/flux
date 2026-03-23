//! flux — composable, type-safe PDE solver framework.
//!
//! flux assembles high-fidelity numerical solvers by composing discrete
//! operators on function spaces defined over simplicial meshes. The central
//! design principle is that incompatible operator compositions (e.g., applying
//! a 1-form operator to a 2-form) are caught at compile time, not runtime.
//!
//! ## Modules
//!
//! - `topology` — simplicial meshes, boundary operators, geometric data
//! - `forms` — discrete k-forms (cochains) parameterized on mesh and degree
//! - `math` — sparse linear algebra primitives (CSR matrices)
//! - `io` — VTK `.vtu` serialization for ParaView visualization

pub const forms = @import("forms/cochain.zig");
pub const io = @import("io/vtk.zig");
pub const math = @import("math/sparse.zig");
pub const operators = struct {
    pub const exterior_derivative = @import("operators/exterior_derivative.zig");
    pub const hodge_star = @import("operators/hodge_star.zig");
};
pub const topology = @import("topology/mesh.zig");

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
