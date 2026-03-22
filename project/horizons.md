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
**Enablers:** None — this can be done now.
**Source:** Ideation 2026-03-22

---

## Dimensionful equation types

**Layer:** architecture
**Constraint on current work:** Keep scalar types parametric (see above). Do not
hardcode `f64` in operator signatures. The eventual design would allow
`Equation(.dimensionful)` vs `Equation(.dimensionless)` so unit checking is
opt-in.
**Enablers:** Generic scalar cochains. Possibly a standalone Zig units library.
**Source:** Ideation 2026-03-22

---

## Pluggable time integrators

**Layer:** architecture
**Constraint on current work:** The M3 leapfrog integrator should not bake
integration logic into the physics module. Separate the staging/stepping
concerns from the operator evaluation. State types should satisfy a comptime
interface (analogous to a Rust trait) that the integrator checks.
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
