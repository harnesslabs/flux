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
