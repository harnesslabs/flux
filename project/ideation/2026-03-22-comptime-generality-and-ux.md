# Ideation: Comptime Generality and Usability

Date: 2026-03-22

## Context

Brainstorming session to identify architectural patterns and design horizons
that flux should be aware of during epoch 1 development. The goal was to
pressure-test raw ideas against the vision and capture validated ones in
`project/horizons.md` so they inform current interface design without becoming
premature commitments.

## Ideas explored

### Intrinsic/extrinsic PDEs and general metrics

**Layer:** vision
**Vision alignment:** extends
**Summary:** The current vision assumes a fixed mesh with flat Euclidean metric.
Intrinsic PDEs (stress/strain, geometric flows) make the mesh a dynamical
variable; general Riemannian and semi-Riemannian metrics make the Hodge star
metric-dependent. The exterior derivative `d` is purely topological and
unaffected — metric dependence is isolated in `★` and inner products.
**Key tension:** Supporting general metrics means `★` is no longer a simple
diagonal matrix. The clean "operators are stateless maps" picture gets muddier
when operators depend on evolving geometry.
**Verdict:** pursue (as design horizon)
**Next action:** Vision amendment added ("Geometric generality as a design
horizon"). Phantom types like `Metric(.flat)` to be introduced when the Hodge
star is implemented (M2). No code for general metrics until a concrete problem
demands it.

### Pluggable time integrators

**Layer:** architecture
**Vision alignment:** aligned
**Summary:** Integrators should be parameterized on method and state type. State
types satisfy a comptime interface (analogous to Rust traits) — leapfrog
requires `staggeredFields()`, RK4 doesn't. The compiler tells you what's
missing.
**Key tension:** Staggered integrators (leapfrog) have a different state shape
than single-snapshot integrators. The interface must accommodate both without
forcing all states to support staggering.
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added with code sketch. First concrete
integrator is M3 leapfrog — design it with the trait interface in mind.

### Dimensionful equation types

**Layer:** architecture
**Vision alignment:** extends
**Summary:** Compile-time dimensional analysis via SI exponents. Opt-in:
`Equation(.dimensionful)` vs `Equation(.dimensionless)`. The dimensionless
case uses raw `f64` with zero overhead. Could be a standalone Zig library.
**Key tension:** Pervasive — if adopted, it touches every scalar in every
operator. Must commit early (generic scalar types) or face a painful retrofit.
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added. Immediate action: make `Cochain` generic
over scalar type (`Cochain(mesh, k, Scalar)` with `Scalar = f64` default).
This also enables `f32` for fast prototyping.

### Operator splitting (Strang and beyond)

**Layer:** feature
**Vision alignment:** aligned
**Summary:** Strang splitting is a meta-integrator that composes two
sub-integrators: `A(dt/2), B(dt), A(dt/2)`. Falls out naturally if the
integrator interface allows callable sub-integrators with arbitrary dt.
**Key tension:** Sub-integrators must not assume they own the full timestep
or the full state.
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added with code sketch. No implementation
until multi-physics use case demands it (likely epoch 2+).

### Domain decomposition (Schwarz + multi-physics)

**Layer:** architecture
**Vision alignment:** aligned (Schwarz), extends (multi-physics coupling)
**Summary:** Schwarz-type methods fall out of mesh views + operator restriction.
Multi-physics coupling (different PDEs on different subdomains) requires
interface operators as first-class nodes. The solver graph may need cycles
(fixed-point coupling), not just DAGs. Mesh partitioning for parallel execution
(MPI) is deferred infrastructure, not a current concern.
**Key tension:** The vision describes the solver graph as a DAG. Multi-physics
coupling may require iteration over a subgraph, making it cyclic.
**Verdict:** pursue (horizon entries for mesh views, Schwarz, multi-physics)
**Next action:** Three horizon entries added. The DAG-vs-cyclic question is
flagged but not resolved — revisit when designing the solver graph abstraction.

### Adaptive methods via comptime strategy selection

**Layer:** constraint
**Vision alignment:** extends
**Summary:** Comptime parameters select adaptation strategy — and with it, the
storage layout, allocator choice, and recomputation policy. `Mesh(.static)` uses
a fixed-buffer allocator; `Mesh(.adaptive)` uses a growth-capable one. The goal:
`Simulation(f32, Mesh(.static), TimeStepper(.fixed), ...)` to
`Simulation(f64, Mesh(.adaptive), TimeStepper(.adaptive), ...)` is a one-line
change.
**Key tension:** Adaptive mesh refinement changes topology, invalidating
cochains, operators, and boundary matrices downstream. Prolongation is
degree-dependent and non-trivial.
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added with concrete constraints on current code
(don't assume fixed topology, don't hardcode operator caching).

### Compositional API and declarative system specification

**Layer:** ux
**Vision alignment:** aligned (compositional), extends (declarative)
**Summary:** Two layers: a compositional API for direct operator chaining
(implementation layer) and a declarative `System` struct for describing PDE
systems at comptime (convenience layer). The declarative form compiles to the
compositional API and supports inspection, consistency checking, and
mathematical display.
**Key tension:** The two layers can diverge if not careful — the declarative
spec must be a strict subset of what the compositional API can express.
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added with code sketches for both layers.
Compositional API emerges naturally from M2 operator work. Declarative spec
is post-M3.

### Structured simulation inputs (JSON + Zig)

**Layer:** ux
**Vision alignment:** extends
**Summary:** JSON manifest for simulation configuration; Zig functions for
boundary conditions and source terms, referenced by name from JSON. No
string-based expression parsers. The JSON is configuration; the Zig code is
behavior.
**Key tension:** The boundary between "configuration" and "behavior" is
fuzzy. Where do material properties go? Initial conditions?
**Verdict:** pursue (horizon entry)
**Next action:** Horizon entry added. No implementation until post-M3.

## Connections

Several of these ideas reinforce each other in a specific way: **comptime
parameters as strategy selectors** is the unifying pattern.

- Generic scalar types enable dimensionful equations
- Pluggable integrators enable operator splitting
- Mesh views enable Schwarz decomposition and multi-physics coupling
- Adaptive strategy parameters compose with all of the above

The composability is multiplicative: each comptime parameter is an independent
axis, and the compiler assembles the right implementation for any combination.
This is the core architectural bet — that Zig's comptime can make this
composition zero-cost and type-safe.

The compositional API and declarative spec sit *above* all of this. They are
the user-facing layer that makes the comptime machinery accessible.

## Open questions

- **DAG vs. cyclic solver graph:** Multi-physics coupling may require
  fixed-point iteration, not just acyclic composition. How does this affect the
  solver graph abstraction?
- **Dimensional analysis: in-tree or dependency?** The units library is
  orthogonal to FEEC/DEC. It could be a standalone Zig package. Decision can
  wait until the need is concrete.
- **Prolongation operators for AMR:** Mapping cochain values onto a refined
  mesh is degree-dependent and non-trivial. This needs a literature review
  before design work begins.
- **Usability validation:** The "domain expert reads the simulation setup"
  test needs to be applied to M3 Maxwell. If the 20-line main.zig isn't
  readable, the API is wrong.

## Artifacts produced

- `project/vision.md`: Added "Geometric generality as a design horizon" and
  "Usability is a correctness property" as design commitments
- `project/horizons.md`: Created with 10 horizon entries, all with code sketches
- `.claude/CLAUDE.md`: Updated to reference horizons.md
- `.claude/commands/{ideate,epoch,tackle,review}.md`: Updated to read/check
  horizons at appropriate points
- `build.zig`: Added `zig build docs` step for API documentation generation
