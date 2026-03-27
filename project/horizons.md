# Design Horizons

Validated architectural ideas that are not yet scheduled. Current work should
not implement these, but should not design in ways that preclude them.

**Check this document before finalizing any interface.**

---

## Generic scalar type for cochains

**Layer:** architecture
**Constraint on current work:** `Cochain` should be parameterized on its scalar
type (`Cochain(mesh, k, Scalar)` with `Scalar = f64` as default) rather than
hardcoding `f64`. This enables `f32` for performance, and later dimensionful
quantities.

```zig
// Current: hardcoded f64
pub fn Cochain(comptime mesh: MeshType, comptime k: u8) type { ... }

// Horizon: generic scalar with default
pub fn Cochain(comptime mesh: MeshType, comptime k: u8, comptime Scalar: type) type {
    // Scalar = f64 for standard use, f32 for fast prototyping,
    // Quantity(.tesla) for dimensionful simulations
    return struct {
        values: []Scalar,
        // ...
    };
}
```

**Enablers:** None — this can be done now.
**Source:** Ideation 2026-03-22

---

## Dimensionful equation types

**Layer:** architecture
**Constraint on current work:** Keep scalar types parametric (see above). Do not
hardcode `f64` in operator signatures. The eventual design would allow
`Equation(.dimensionful)` vs `Equation(.dimensionless)` so unit checking is
opt-in.

```zig
// Compile-time dimensional analysis via SI exponents (M, L, T, I, ...)
pub fn Quantity(comptime dim: Dimensions) type {
    return struct {
        value: f64,
        // arithmetic checks dimensional compatibility at comptime
        pub fn mul(self: @This(), other: anytype) Quantity(dim.add(@TypeOf(other).dim)) { ... }
    };
}

const Length = Quantity(.{ .length = 1 });
const Time = Quantity(.{ .time = 1 });
const Velocity = Quantity(.{ .length = 1, .time = -1 });

// Opt-in: dimensionless just uses raw f64
const maxwell_fast = System(.dimensionless, ...); // scalars are f64
const maxwell_safe = System(.dimensionful, ...);  // scalars are Quantity(...)
```

**Enablers:** Generic scalar cochains. Possibly a standalone Zig units library.
**Source:** Ideation 2026-03-22

---

## Pluggable time integrators

**Layer:** architecture
**Constraint on current work:** The M3 leapfrog integrator should not bake
integration logic into the physics module. Separate the staging/stepping
concerns from the operator evaluation. State types should satisfy a comptime
interface (analogous to a Rust trait) that the integrator checks.

```zig
// The integrator is parameterized on method and state type
pub fn Integrator(comptime method: IntegrationMethod, comptime State: type) type {
    // Verify State satisfies the required interface at comptime
    comptime {
        if (!@hasDecl(State, "evaluateOperator"))
            @compileError("State must declare evaluateOperator");
        if (method == .leapfrog and !@hasDecl(State, "staggeredFields"))
            @compileError("Leapfrog requires State to declare staggeredFields");
    }
    return struct {
        pub fn step(state: *State, dt: f64) StepResult(method) { ... }
    };
}

// Users build their state struct, declare the required functions
const MaxwellState = struct {
    E: Cochain(mesh, 1, f64),
    B: Cochain(mesh, 2, f64),
    pub fn evaluateOperator(self: *@This(), ...) void { ... }
    pub fn staggeredFields(self: *@This()) ... { ... }
};

const stepper = Integrator(.leapfrog, MaxwellState);
```

**Enablers:** M2 operators, M3 Maxwell simulation (first concrete integrator).
**Source:** Ideation 2026-03-22

---

## Operator splitting (Strang and beyond)

**Layer:** feature
**Constraint on current work:** Integrators must be composable — a
meta-integrator like Strang splitting takes two sub-integrators and interleaves
their steps (`A(dt/2), B(dt), A(dt/2)`). This means the integrator interface
cannot assume it owns the full timestep or the full state. Sub-integrators must
be callable with arbitrary dt and must not hardcode step-size logic.

```zig
// Strang splitting is a meta-integrator: it composes two sub-integrators
pub fn StrangSplit(comptime IntA: type, comptime IntB: type) type {
    return struct {
        a: IntA,
        b: IntB,
        pub fn step(self: *@This(), state: *IntA.State, dt: f64) void {
            self.a.step(state, dt / 2.0);
            self.b.step(state, dt);
            self.a.step(state, dt / 2.0);
        }
    };
}

// Usage: split advection from diffusion
const split = StrangSplit(
    Integrator(.explicit, AdvectionOp),
    Integrator(.implicit, DiffusionOp),
);
```

**Enablers:** Pluggable time integrators with a generic interface. First real
use case is likely multi-physics (advection/diffusion split, field/particle
split).
**Source:** Ideation 2026-03-22

---

## Mesh views and submesh references

**Layer:** architecture
**Constraint on current work:** Mesh data structures should not assume they own
all their data exclusively. Leave room for a view/slice type that references a
parent mesh's storage for a subset of entities.
**Enablers:** Stable mesh data layout.
**Source:** Ideation 2026-03-22

---

## Schwarz-type domain decomposition

**Layer:** architecture
**Constraint on current work:** None immediately — this is a solver-graph
pattern that falls out of mesh views + operator restriction.
**Enablers:** Mesh views, operator restriction to submeshes.
**Source:** Ideation 2026-03-22

---

## Multi-physics coupling at interfaces

**Layer:** architecture
**Constraint on current work:** The solver graph model (currently described as a
DAG in `vision.md`) may need to support cycles for fixed-point coupling between
subdomains with different physics. Keep this in mind when designing the solver
graph abstraction — don't assume acyclicity structurally.
**Enablers:** Interface operators (trace, projection, interpolation) as
first-class nodes. Mesh views.
**Source:** Ideation 2026-03-22

---

## Adaptive methods via comptime strategy selection

**Layer:** constraint
**Constraint on current work:** Design meshes, integrators, and operators so that
a comptime parameter selects the adaptation strategy — and with it, the storage
layout, allocator choice, and recomputation policy.

Concrete constraints on current code:
- **Meshes:** Do not assume topology is fixed forever. Don't store operator
  matrices as persistent state that is never rebuilt. Use entity indices rather
  than raw array positions where cost is acceptable.
- **Integrators:** Step functions should return what happened (actual dt, error
  estimate, accepted/rejected), not just the new state. A `TimeStepper(.fixed)`
  returns only state; `TimeStepper(.adaptive)` returns the full diagnostic. The
  comptime parameter determines which fields exist.
- **Cochains:** Do not assume cochain length equals mesh entity count for all
  time. Adaptive mesh refinement requires prolongation (mapping cochain values
  onto a refined mesh), which is degree-dependent.
- **Allocators:** `Mesh(.static)` can use a fixed-buffer allocator (one
  allocation, cache-perfect layout). `Mesh(.adaptive)` needs a growth-capable
  allocator. The comptime mesh parameter should influence which allocator
  strategy is appropriate — Zig's explicit allocators make this a natural
  extension rather than a fight against defaults.

The goal: changing `Simulation(f32, Mesh(.static), TimeStepper(.fixed), ...)` to
`Simulation(f64, Mesh(.adaptive), TimeStepper(.adaptive), ...)` is a one-line
change that selects a coherent, optimized implementation strategy. The user
declares intent; the compiler assembles the machinery.

**Enablers:** Generic scalar cochains, pluggable integrators, stable mesh
abstraction.
**Source:** Ideation 2026-03-22

---

## Compositional API and declarative system specification

**Layer:** ux
**Constraint on current work:** Operators (`d`, `★`, compositions) should be
designed as composable units from day one — this is already in the vision. The
further horizon is a declarative `System` struct that describes a PDE system at
comptime (fields, equations, coupling) and compiles down to the compositional
API. The declarative form enables inspection, consistency checking, and
formatted display of the system in mathematical notation.

Well-formedness invariants (all comptime-checkable):
- Every field referenced in an equation exists in the system
- Every operator is dimensionally compatible with its input
- The system is not underdetermined

The compositional API is the implementation layer; the declarative spec is a
convenience that compiles down to it. Power users compose directly; standard
setups use the declarative form.

```zig
// Compositional API (implementation layer) — direct operator chaining
const curl_E = d(1).apply(E);
const dB_dt = curl_E.negate();
const J_curl_H = star_inv(1).compose(d(0)).compose(star(2)).apply(B);

// Declarative spec (convenience layer) — compiles to the above at comptime
const maxwell = System{
    .fields = .{
        Field("E", .{ .mesh = mesh, .degree = 1, .duality = .primal }),
        Field("B", .{ .mesh = mesh, .degree = 2, .duality = .primal }),
    },
    .equations = .{
        .{ .ddt = "B", .rhs = .{ .negate, .{ .d, "E" } } },
        .{ .ddt = "E", .rhs = .{ .sub, .{ .star_inv, .{ .d, .{ .star, "B" } } }, "J" } },
    },
};

// Declarative form supports inspection and display
maxwell.printSystem(); // prints: ∂B/∂t = −dE,  ∂E/∂t = ★⁻¹d★B − J
```

**Enablers:** Typed cochains, operator composition (M2), at least one working
simulation (M3) to validate the API shape against real physics code.
**Source:** Ideation 2026-03-22

---

## Structured simulation inputs (JSON + Zig)

**Layer:** ux
**Constraint on current work:** None immediately — this is a frontend concern
that sits above the operator/solver layer. But do not bake parameter handling
into simulation code in ways that would be hard to externalize later.

Design intent: a JSON manifest describes the simulation configuration (mesh
file path, timestep, output settings, named boundary conditions and sources).
Boundary conditions and source terms that require expressions are Zig functions
registered by name — the JSON references them, not inlines them. The JSON is
configuration; the Zig code is behavior. No string-based expression parsers.

```json
{
  "mesh": "cavity.obj",
  "time": { "dt": 1e-3, "steps": 1000 },
  "boundary_conditions": { "walls": { "type": "pec" } },
  "sources": { "dipole": { "function": "dipole_source" } }
}
```

**Enablers:** A working simulation (M3), stable public API.
**Source:** Ideation 2026-03-22

---

## Native quad/polygon cell support

**Layer:** architecture
**Constraint on current work:** The `.obj` reader must produce a `RawMesh`
intermediate representation that preserves polygon faces of any valence
(triangles, quads, n-gons) *before* an explicit, separable `triangulate()` pass
converts to simplicial form. The parser is format-aware; the triangulator is
topology-aware. They must remain separate passes so that future work can consume
polygon cells directly without reworking the I/O layer.

The deeper constraint: all Whitney-based operators (convergent Hodge star ★₁,
field reconstruction) assume simplicial complexes. Supporting quads natively
requires a parallel code path with tensor-product basis functions (FEM-style Q1
elements) — a different discretization theory. The topological operators (∂, d)
are already cell-shape-agnostic and would work on quads immediately.

```zig
// Current (M2): RawMesh → triangulate() → simplicial Mesh → operators
// Horizon:      RawMesh → Mesh(CellType.quad) → quad-aware operators
//               RawMesh → Mesh(CellType.simplex) → simplicial operators (existing)

const CellType = enum { simplex, quad, polygon };

// Mesh parameterized on cell type — operators dispatch accordingly
pub fn Mesh(comptime n: u8, comptime dim: u8, comptime cell: CellType) type { ... }
```

**What works on quads today (topological, no simplicial assumption):**
- Boundary operators `∂ₖ` — a quad has 4 boundary edges
- Exterior derivative `d` — just `∂ᵀ`, topology only
- Diagonal Hodge star — needs per-cell volume formulas, but quad area / dual
  edge length formulas are straightforward

**What requires new theory for quads:**
- Whitney/Galerkin mass matrix — Whitney forms use barycentric coordinates on
  simplices. Quads need bilinear (Q1) basis functions.
- Any operator that interpolates within cells (field reconstruction, evaluation)

**Enablers:** `RawMesh` intermediate representation (M2 #86), dimension-generic
topology (M2 #80), stable operator trait interfaces (M1 #75).
**Source:** M2 milestone planning 2026-03-26
