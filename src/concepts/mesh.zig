//! Comptime mesh concept.
//!
//! `MeshConcept` validates that a type provides the minimal interface that
//! DEC operators need to function: entity counts, boundary operators, and
//! topological/embedding dimension metadata.
//!
//! This decouples operators from the concrete `topology.Mesh(n)` type,
//! enabling dimension-generic operators and alternative mesh implementations
//! (e.g., mesh views, imported meshes) without code duplication.
//!
//! Conforming types must declare:
//!   - `pub const dimension: usize`              — embedding dimension
//!   - `pub const topological_dimension: usize`   — mesh dimension (2 for surfaces)
//!   - `pub fn num_vertices(self) u32`
//!   - `pub fn num_edges(self) u32`
//!   - `pub fn num_faces(self) u32`
//!   - `pub fn boundary(self, comptime k) BoundaryMatrix` — boundary operator ∂ₖ

const std = @import("std");
const testing = std.testing;

// ═══════════════════════════════════════════════════════════════════════════
// MeshConcept — comptime validation
// ═══════════════════════════════════════════════════════════════════════════

/// Validate that `M` satisfies the Mesh concept at compile time.
///
/// A conforming mesh type must declare:
///   - `pub const dimension: usize`             — embedding dimension ℝⁿ
///   - `pub const topological_dimension: usize`  — intrinsic dimension (2 for surfaces)
///   - `pub fn num_vertices(self) u32`
///   - `pub fn num_edges(self) u32`
///   - `pub fn num_faces(self) u32`
///   - `pub fn boundary(self, comptime k: comptime_int) BoundaryMatrix`
///
/// Produces a descriptive `@compileError` on violation.
pub fn MeshConcept(comptime M: type) void {
    // 1. M must declare an embedding dimension.
    if (!@hasDecl(M, "dimension")) {
        @compileError("MeshConcept requires a 'pub const dimension' declaration — " ++
            "the embedding dimension (e.g., 2 for ℝ²)");
    }

    // 2. M must declare a topological dimension.
    if (!@hasDecl(M, "topological_dimension")) {
        @compileError("MeshConcept requires a 'pub const topological_dimension' declaration — " ++
            "the intrinsic mesh dimension (e.g., 2 for surfaces)");
    }

    // 3. Entity count accessors: num_vertices, num_edges, num_faces.
    // Each must be a method taking self and returning u32.
    inline for (.{ "num_vertices", "num_edges", "num_faces" }) |name| {
        if (!@hasDecl(M, name)) {
            @compileError("MeshConcept requires a 'pub fn " ++ name ++ "(self) u32' declaration");
        }
        const decl_info = @typeInfo(@TypeOf(@field(M, name)));
        if (decl_info != .@"fn") {
            @compileError("MeshConcept: '" ++ name ++ "' must be a function");
        }
        const fn_info = decl_info.@"fn";
        if (fn_info.params.len != 1) {
            @compileError("MeshConcept: '" ++ name ++ "' must take exactly 1 parameter (self)");
        }
        const ret = fn_info.return_type orelse
            @compileError("MeshConcept: '" ++ name ++ "' must have a known return type");
        if (ret != u32) {
            @compileError("MeshConcept: '" ++ name ++ "' must return u32");
        }
    }

    // 4. M must declare a boundary function.
    if (!@hasDecl(M, "boundary")) {
        @compileError("MeshConcept requires a 'pub fn boundary(self, comptime k: comptime_int)' declaration — " ++
            "returns the boundary operator ∂ₖ");
    }
    const boundary_info = @typeInfo(@TypeOf(M.boundary));
    if (boundary_info != .@"fn") {
        @compileError("MeshConcept: 'boundary' must be a function");
    }
    // boundary takes (self, comptime k) — 2 params.
    if (boundary_info.@"fn".params.len != 2) {
        @compileError("MeshConcept: 'boundary' must take exactly 2 parameters (self, comptime k)");
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests — MeshConcept
// ═══════════════════════════════════════════════════════════════════════════

/// A minimal conforming mock mesh for testing the concept.
const MockMesh = struct {
    pub const dimension = 2;
    pub const topological_dimension = 2;

    vertex_count: u32,
    edge_count: u32,
    face_count: u32,

    pub fn num_vertices(self: @This()) u32 {
        return self.vertex_count;
    }

    pub fn num_edges(self: @This()) u32 {
        return self.edge_count;
    }

    pub fn num_faces(self: @This()) u32 {
        return self.face_count;
    }

    /// Stub boundary operator — concept only checks the declaration exists.
    pub fn boundary(self: @This(), comptime k: comptime_int) MockBoundaryMatrix {
        _ = self;
        _ = k;
        unreachable;
    }
};

const MockBoundaryMatrix = struct {
    n_rows: u32,
    n_cols: u32,
};

test "MeshConcept accepts a conforming type" {
    comptime MeshConcept(MockMesh);
}

test "MeshConcept accepts a type with extra declarations" {
    const ExtendedMesh = struct {
        pub const dimension = 3;
        pub const topological_dimension = 2;

        pub fn num_vertices(_: @This()) u32 {
            return 0;
        }
        pub fn num_edges(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
        pub fn boundary(_: @This(), comptime k: comptime_int) MockBoundaryMatrix {
            _ = k;
            unreachable;
        }

        /// Extra — concept should not reject.
        pub fn volume(_: @This()) f64 {
            return 0;
        }
    };
    comptime MeshConcept(ExtendedMesh);
}

test "MeshConcept accepts the real Mesh(2) type" {
    const topology = @import("../topology/mesh.zig");
    comptime MeshConcept(topology.Mesh(2));
}

// ── Negative tests (compile-time rejection) ──────────────────────────────

test "MeshConcept rejects type missing dimension" {
    const NoDim = struct {
        pub const topological_dimension = 2;
        pub fn num_vertices(_: @This()) u32 {
            return 0;
        }
        pub fn num_edges(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
        pub fn boundary(_: @This(), comptime k: comptime_int) MockBoundaryMatrix {
            _ = k;
            unreachable;
        }
    };

    // comptime MeshConcept(NoDim);
    // expected: @compileError("MeshConcept requires a 'pub const dimension' declaration ...")
    _ = NoDim;
}

test "MeshConcept rejects type missing topological_dimension" {
    const NoTopoDim = struct {
        pub const dimension = 2;
        pub fn num_vertices(_: @This()) u32 {
            return 0;
        }
        pub fn num_edges(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
        pub fn boundary(_: @This(), comptime k: comptime_int) MockBoundaryMatrix {
            _ = k;
            unreachable;
        }
    };

    // comptime MeshConcept(NoTopoDim);
    // expected: @compileError("MeshConcept requires a 'pub const topological_dimension' declaration ...")
    _ = NoTopoDim;
}

test "MeshConcept rejects type missing num_vertices" {
    const NoVerts = struct {
        pub const dimension = 2;
        pub const topological_dimension = 2;
        pub fn num_edges(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
        pub fn boundary(_: @This(), comptime k: comptime_int) MockBoundaryMatrix {
            _ = k;
            unreachable;
        }
    };

    // comptime MeshConcept(NoVerts);
    // expected: @compileError("MeshConcept requires a 'pub fn num_vertices(self) u32' declaration ...")
    _ = NoVerts;
}

test "MeshConcept rejects type missing num_edges" {
    const NoEdges = struct {
        pub const dimension = 2;
        pub const topological_dimension = 2;
        pub fn num_vertices(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
        pub fn boundary(_: @This(), comptime k: comptime_int) MockBoundaryMatrix {
            _ = k;
            unreachable;
        }
    };

    // comptime MeshConcept(NoEdges);
    // expected: @compileError("MeshConcept requires a 'pub fn num_edges(self) u32' declaration ...")
    _ = NoEdges;
}

test "MeshConcept rejects type missing boundary" {
    const NoBoundary = struct {
        pub const dimension = 2;
        pub const topological_dimension = 2;
        pub fn num_vertices(_: @This()) u32 {
            return 0;
        }
        pub fn num_edges(_: @This()) u32 {
            return 0;
        }
        pub fn num_faces(_: @This()) u32 {
            return 0;
        }
    };

    // comptime MeshConcept(NoBoundary);
    // expected: @compileError("MeshConcept requires a 'pub fn boundary(...)' declaration ...")
    _ = NoBoundary;
}
