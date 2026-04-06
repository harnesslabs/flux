# Epoch 2 Decision Log

<!-- Append entries with /decide. Format: ## YYYY-MM-DD: <title> -->

## 2026-03-29: TimeStepStrategy concept uses structural signature checking

**Decision:** `TimeStepStrategy(S)` validates `S.step` via `@typeInfo` — checking parameter
count, parameter types, and that the return is an error union with `void` payload — rather
than comparing `@TypeOf(S.step) == fn(Allocator, *State, f64) anyerror!void`.

**Alternatives considered:**
1. Exact type equality (`@TypeOf(S.step) == expected`): rejected because Zig distinguishes
   inferred error sets from `anyerror`. A function declared `fn(...) !void` gets a concrete
   error set, not `anyerror`, so no real implementation would match the exact type.
2. Duck typing (no concept check, just `anytype`): rejected because silent failures at
   call sites are harder to debug than a concept-level compile error.

**Rationale:** Structural checking accepts any error set in the return type, which is the
correct generality — an integrator's error set depends on which operators it calls. The
concept guarantees the parameter contract (allocator, state pointer, dt) while staying
permissive on the error set.

## 2026-03-29: Two-layer TimeStepper design — concept + generic wrapper

**Decision:** Split into `TimeStepStrategy` (comptime concept) and `TimeStepper(Strategy)`
(generic wrapper). Users define a strategy struct with the raw physics/numerics, then
pass it to `TimeStepper(...)` which validates the contract and provides a uniform
integrator with step counting.

**Alternatives considered:**
1. Single `TimeStepper` concept with no wrapper: simpler, but forces every consumer
   (runners, composition) to re-implement step counting and any shared integrator logic.
2. Inheritance-style base struct: not idiomatic Zig. Comptime generics with concept
   validation is the Zig equivalent of Rust's `impl Trait` pattern.

**Rationale:** The two-layer design separates concerns — strategy authors define physics,
the wrapper provides infrastructure (step counting now, diagnostics/observers later).
This matches the "pluggable time integrators" horizon and enables future composition
(Strang splitting wraps two `TimeStepper` instances).

## 2026-03-29: OperatorConcept deferred — no consumer yet

**Decision:** Do not introduce `OperatorConcept` now. Existing operators are free
functions using `anytype`; Zig's structural type system already gives compile errors
on degree/duality mismatches via the cochain type. `compose.chain` threads return
types at comptime without needing explicit `InputType`/`OutputType` metadata.

**What was tried and removed:**
An `OperatorConcept(Op)` validator was prototyped along with wrapper types
(`ExteriorDerivativeOperator`, `HodgeStarOperator`, etc.) that packaged free
functions into concept-conforming structs. The wrappers were thin delegates with
no consumer — no generic code dispatched on `OperatorConcept`-conforming types.

**When to revisit:**
When there is actual generic infrastructure that needs to inspect operator
input/output types at comptime — e.g., declarative system specs (the `System`
horizon), solver graphs, or operator algebra composition that goes beyond
`compose.chain`. Design the concept to fit the consumer, not the other way around.

**Rationale:** A concept without a consumer is premature abstraction. `MeshConcept`
has a clear near-term consumer (mesh views, dimension-generic meshes). `TimeStepper`
has a consumer (the `Runner`). Operators don't have one yet — the `anytype` + cochain
type system already provides the guarantees users need for interoperability.

## 2026-03-29: MeshConcept — minimal interface, geometric data excluded

**Decision:** `MeshConcept(M)` requires only: `dimension`, `topological_dimension`,
`num_vertices/edges/faces(self) u32`, and `boundary(self, comptime k)`. Geometric data
(dual areas, edge lengths, Whitney mass matrix) is *not* part of the concept.

**Alternatives considered:**
1. Include geometric accessors (dual_area, edge_length, etc.): rejected because these
   are operator-specific. The Hodge star needs dual volumes; the exterior derivative
   doesn't. Requiring all geometric data in the concept would over-constrain alternative
   mesh implementations (mesh views, imported meshes) that may not have all data computed.
2. Separate `GeometricMesh` concept extending `MeshConcept`: considered for the future
   but premature now. When M2's dimension-generic operators land, the needed geometric
   interface will be clearer.

**Rationale:** The minimal concept captures what *all* DEC operators need: topology (entity
counts, boundary operators) and dimension metadata. Geometric data varies by operator and
should be checked by the operators themselves, not the mesh concept.

## 2026-03-29: Concept naming convention — `*Concept` suffix

**Decision:** Comptime concept validators use the `*Concept` suffix: `MeshConcept`.
This is slightly inconsistent with the existing `TimeStepStrategy` name (which serves
the same role) but more descriptive.

**Alternatives considered:**
1. `*Strategy` suffix for all: rejected because "strategy" implies a behavioral pattern
   (Strategy Pattern), not a type contract. "Concept" is the standard C++ term for
   compile-time type constraints and maps well to Zig's comptime validation.
2. Rename `TimeStepStrategy` → `TimeStepConcept` now: deferred to avoid scope creep.
   Should be done in a follow-up issue to unify naming.

**Rationale:** `*Concept` is the clearest name for "comptime function that validates a
type satisfies an interface." Unifying `TimeStepStrategy` is tracked as a follow-up.

## 2026-03-29: Benchmark baselines generated on CI, not committed locally

**Decision:** `bench/baselines.json` is generated and committed by a dedicated
`update-baselines.yml` GitHub Actions workflow running on `ubuntu-latest`. Local
`zig build bench` runs without `--check` by default (informational only). The CI
`bench` job runs `--check` to detect regressions against the CI-generated baselines.

**Alternatives considered:**
1. Commit local (developer machine) baselines: rejected because absolute timings
   vary wildly across machine architectures (ARM vs x86, different clock speeds).
   A 20% threshold would produce false positives on every machine transition.
2. Ratio-based baselines (normalize to a reference benchmark): considered but adds
   complexity. Ratios are fragile if the reference benchmark's absolute time changes
   due to system load rather than code changes.
3. No committed baselines — benchmarks are purely informational: rejected because
   the M1 acceptance criterion requires "bounded CI checks."

**Rationale:** Same-machine-class comparison (ubuntu-latest → ubuntu-latest) keeps
the 20% regression threshold meaningful. The `update-baselines` workflow runs
automatically when operator code changes on main, keeping baselines fresh.

## 2026-03-30: Generic integrators live in `src/integrators/`, not `src/operators/`

**Decision:** Time integrators (leapfrog, forward Euler, future backward Euler)
live in `src/integrators/`. The issue suggested `src/operators/` as an alternative.

**Alternatives considered:**
1. `src/operators/` — rejected because `operators/` contains *spatial* DEC operators
   (d, ★, Δ). Time integrators are a fundamentally different concern — they compose
   spatial operators into time-stepping patterns. Mixing them blurs the boundary.
2. `src/time_stepper.zig` (expand the existing file) — rejected because the existing
   file defines the `TimeStepStrategy` concept and `TimeStepper` wrapper. Concrete
   integrators are consumers of that concept, not part of its definition.

**Rationale:** Separate directory keeps the spatial/temporal distinction clear and
gives a natural home for future integrators (backward Euler for #93 heat equation,
Strang splitting from horizons).

## 2026-04-01: Assembled Laplacian is a standalone operator object, not mesh-owned cache

**Decision:** Represent the pre-assembled scalar Laplacian as an explicit
operator value in `src/operators/laplacian.zig` with its own lifetime:
`assemble_primal_0_laplacian(allocator, mesh) -> Primal0Laplacian(MeshType)`.
The object owns the stiffness matrix `S = D₀ᵀ M₁ D₀` and the diagonal `★₀⁻¹`
data needed for application.

**Alternatives considered:**
1. Add a Laplacian cache field directly to `Mesh`: rejected because it pushes
   solver/operator cache policy into topology ownership, making `Mesh` harder
   to evolve toward mesh views and alternate operator families.
2. Hide assembly inside `laplacian(...)` and rebuild every call: rejected
   because it does not satisfy the issue's repeated-application goal and keeps
   allocation cost on the hot path.

**Rationale:** A standalone operator object makes ownership explicit, keeps the
topology layer focused on mesh data, and leaves room for multiple assembled
operators over the same mesh (scalar Laplacian, vector Laplacian, metric
variants, solver preconditioners) without mutating mesh state. This is closer
to the compositional-operator direction in `horizons.md` than embedding caches
inside the mesh itself.

## 2026-04-02: Assembled operators live in `OperatorContext`, not as standalone public wrappers

**Decision:** Introduce `OperatorContext(MeshType)` as the owner of assembled
DEC operators for one mesh and one problem. The context exposes
`withLaplacian(k)` and `laplacian(k)`; the previously-added public Laplacian
convenience wrappers are removed instead of kept in parallel.

**Alternatives considered:**
1. Keep standalone public wrappers alongside the context: rejected because it
   creates two ways to express the same operation, and the one-shot path hides
   assembly cost in exactly the place hot-loop code should be explicit.
2. Move caches directly onto `Mesh`: rejected again because topology/geometry
   ownership and problem-specific operator assembly are different concerns.

**Rationale:** The mesh owns reusable geometric facts; the operator context owns
only the assembled state actually required by the current system. This keeps the
API honest about setup cost, avoids parallel interfaces in a pre-release codebase,
and gives the future PDE/system builder a natural home: the builder requests the
operators it needs from the context instead of mutating the mesh or relying on
hidden convenience assembly.

## 2026-04-03: SIMD cochain arithmetic uses private slice kernels, not new public methods

**Decision:** Keep the caller-facing cochain API unchanged:
`Cochain.add`, `.scale`, `.negate`, and `.inner_product` remain the only
arithmetic entry points. Implement SIMD through private slice kernels inside
`src/forms/cochain.zig`, with scalar reference kernels kept alongside them for
tests and benchmark baselines.

**Alternatives considered:**
1. Add explicit public SIMD methods (`addSimd`, `innerProductSimd`, etc.):
   rejected because it creates parallel APIs for the same algebraic operations,
   forcing callers to care about an implementation detail.
2. Keep only one SIMD implementation with no scalar reference kernel:
   rejected because the issue requires measurable benchmark comparison against a
   scalar baseline, and the tail-handling path is easier to validate against a
   trusted scalar implementation than against hand-derived expected arrays.

**Rationale:** This keeps one obvious way to express cochain arithmetic while
still giving the module an internal performance seam. The decision also stays
compatible with the scalar-type horizon: when `Cochain(..., Scalar)` lands, the
private kernels can be generalized or specialized without freezing a public
SIMD-specific surface today.

## 2026-04-03: Incidence compression prototypes a sign-bit CSR wrapper before any mesh integration

**Decision:** Prototype compressed incidence storage as a dedicated
`PackedIncidenceMatrix` plus `PackedIncidenceSigns` in `src/math/sparse.zig`,
benchmark it against `CsrMatrix(i8)`, and defer replacing mesh boundary
operators until the benchmark shows a real SpMV win.

**Alternatives considered:**
1. Integrate a 2-bit packed ternary encoding directly into `CsrMatrix(i8)`:
   rejected for the first pass because CSR already omits zeros, so storing
   `{−1,+1}` as one sign bit per nonzero is denser, and changing `row().vals`
   semantics up front would force avoidable API churn across operators.
2. Replace `BoundaryMatrix` immediately with a compressed type: rejected until
   the packed path proves it helps the hot loop. The issue explicitly allows
   not shipping the optimization if decode overhead erases the cache benefit.

**Rationale:** The real question in issue #48 is empirical: does smaller
incidence storage make boundary-operator application faster on large meshes?
Separating the storage experiment from mesh/operator integration keeps the API
honest, preserves the current topology/operator boundary, and answers the
performance question before committing the rest of the codebase to a more
involved matrix interface.

## 2026-04-04: Uniform tetrahedral grids use the 6-tet Freudenthal/Kuhn decomposition

**Decision:** `uniform_tetrahedral_grid(...)` decomposes each axis-aligned cube
into six tetrahedra using the Freudenthal/Kuhn simplex decomposition induced by
the coordinate ordering `x ≤ y ≤ z` on the unit cube, then maps that pattern
affinely onto every hexahedral cell in the grid.

**Alternatives considered:**
1. Five-tetrahedron cube split with parity alternation: rejected because the
   global conformity rules are harder to reason about and would force
   checkerboard bookkeeping into what should be a mechanically generated test
   bed.
2. Ad hoc per-cell triangulations chosen for symmetry: rejected because local
   choices make face sharing and orientation consistency across cell boundaries
   harder to audit.

**Rationale:** The 6-tet decomposition gives one obvious global rule: every
cell uses the same simplex pattern derived from the cube's total vertex order.
That makes boundary assembly easier to verify, guarantees conforming shared
faces across adjacent cells, and provides a clean foundation for M2 invariant
tests. The tradeoff is a larger tetrahedron count than a 5-tet scheme, which
is acceptable here because correctness and orientation regularity matter more
than minimizing element count in the initial 3D test-bed constructor.

## 2026-04-04: Interior-degree Hodge stars are stored as a degree-indexed Whitney family on the mesh

**Decision:** Replace the mesh's single `whitney_mass_1`/`preconditioner_1`
fields with a degree-indexed `whitney_operators` family covering every
interior degree `1..n-1`. The Hodge star dispatch remains simple:
boundary degrees `k = 0` and `k = n` use diagonal geometric ratios, while
every interior degree uses the Whitney/Galerkin mass matrix `M_k` for forward
application and PCG for `★⁻¹`, initialized from the diagonal dual/primal ratio.

## 2026-04-05: Observers live beside operators and evaluate explicitly, not inside `TimeStepper`

**Decision:** Add the M3 diagnostic framework as standalone observer values in
`src/operators/observers.zig`, with explicit `evaluate(allocator, state, step)`
methods plus a tuple-based `evaluateAll(...)` helper. Do not extend
`TimeStepper` yet.

**Alternatives considered:**
1. Add observer attachment directly to `TimeStepper`: rejected because there is
   no shared runner API yet, so baking storage and output policy into the
   wrapper would force an early interface around an unknown consumer.
2. Hide allocation inside observer objects so `evaluate(state, step)` stays
   allocator-free: rejected because diagnostics like divergence norm and
   helicity legitimately allocate temporary cochains today. Hiding that would
   violate the project's explicit-allocation rule.

**Rationale:** The current milestone needs reusable diagnostics more urgently
than it needs a universal simulation runner. Standalone observers give every
example and future runner one obvious hook to call after each step, keep
allocation explicit, and preserve the option to let a later runner own the
policy for buffering, formatting, and persistence without first unwinding an
over-eager `TimeStepper` interface.

**Alternatives considered:**
1. Keep a special `whitney_mass_1` path and add a separate 3D-only `k = 2`
   implementation beside it: rejected because it hard-codes today's dimensions
   into the API and duplicates the same ownership/lifetime pattern under
   different names.
2. Keep the old diagonal 3D interior Hodge formulas as public execution paths:
   rejected because on barycentric duals they are preconditioners, not the
   correct Galerkin operator, so preserving them as peers would repeat the
   same mistake issue #70 fixed for 2D `★₁`.

**Rationale:** The geometry/topology layer already owns the data every Whitney
Hodge star needs: primal measures, dual measures, simplex incidence, and
barycenters. Storing interior operators by degree keeps the public surface at
one obvious rule, removes dimension-specific API bloat, and leaves the correct
mathematical specialization in the implementation rather than in caller-visible
field names. Using the diagonal ratio as the PCG preconditioner and initial
guess preserves the DEC intuition as a fast spectral approximation without
pretending it is the operator itself.

## 2026-04-05: Wedge product is defined by Whitney interpolation and de Rham projection

**Decision:** Define the primal-primal wedge as the lowest-order FEEC/Whitney
product induced by interpolation and projection:
`wedge(α, β) := R(Wα ∧ Wβ)`, where `W` is the Whitney map and `R` is the de Rham
map back to simplicial cochains. The implementation computes the equivalent
local simplex formula directly instead of constructing the interpolated forms at
runtime.

**Alternatives considered:**
1. Treat wedge as a purely combinatorial cochain product with no FEEC/Whitney
   interpretation: rejected because flux has already committed to Whitney-based
   operators for correctness on the Hodge-star path, so a nonlinear product with
   a separate mathematical story would fragment the operator stack.
2. Delay wedge until a full higher-order FEEC field layer exists: rejected
   because M3 needs a usable wedge now for vorticity and helicity work, and the
   lowest-order Whitney-induced product is a coherent stepping stone.
3. Require strict associativity on arbitrary cochains: rejected because the
   projected lowest-order Whitney product does not satisfy it in general. The
   correct guarantee is graded commutativity and Leibniz exactly, with
   associativity only for closed forms.

**Rationale:** This choice points the codebase toward the FEEC/Whitney future
instead of enshrining cochains as the only semantic layer, while still exposing
the simple operator surface the current library needs. It also leaves a clean
upgrade path: future work can make `W` and `R` explicit, add richer finite
element form spaces, and represent nonlinear products without immediately
collapsing them back to the lowest-order cochain complex.

## 2026-04-05: Dirichlet Poisson solve uses symmetric row/column elimination on the stiffness system

**Decision:** Enforce Dirichlet boundary conditions for the scalar Poisson solve
by reducing the symmetric stiffness system `S u = ★₀ f` to interior vertices:
remove boundary rows and columns, move the eliminated boundary contribution into
the right-hand side, and solve the reduced SPD system with CG.

**Alternatives considered:**
1. Penalty method: rejected because it introduces a large artificial parameter,
   perturbs the operator being verified, and weakens the convergence claim.
2. Row overwrite only (`u_i = g_i` on boundary rows): rejected because it breaks
   symmetry unless columns are also handled carefully, which is unnecessary
   complexity when we already have a clean reduced formulation.
3. Solve the full Laplacian `Δ₀ u = f` directly: rejected because the assembled
   operator is stored as `Δ₀ = ★₀⁻¹ S`; the stiffness form is the natural SPD
   object for Dirichlet elimination and CG.

**Rationale:** The reduced system preserves the exact discrete operator on the
interior, keeps the matrix SPD for CG, and makes the boundary treatment explicit.
This fits the operator-context architecture: topology owns geometric facts,
operators own assembled algebra, and boundary handling is a solver-layer concern
rather than a mutation of the mesh or Laplacian operator itself.

## 2026-04-05: CG preconditioners use a comptime concept and typed pointer argument

**Decision:** The conjugate-gradient solver accepts the preconditioner as a
typed pointer argument (`null`, `?*const P`, or `*const P`) validated by a
`PreconditionerConcept(P)` comptime check. A conforming preconditioner exposes
`pub fn apply(self: *const P, z: []f64) void`, mutating only the work vector.

**Alternatives considered:**
1. Keep the erased interface (`fn ([]f64, *const anyopaque) void` plus a
   separate context pointer): rejected because it discards the preconditioner
   type exactly where Zig can validate it, and every call site must keep the
   function and context in sync manually.
2. Add a second typed solve entry point but keep the erased `solve(...)`
   public as well: rejected because this is pre-release code with no real
   compatibility pressure. Two public solver entry points for one algorithm
   would violate the "one obvious way" rule.

**Rationale:** The typed-pointer design keeps solver call sites honest, matches
the existing concept-driven style already used for time stepping and meshes,
and makes Jacobi preconditioning a normal Zig value instead of a function/context
pair. Requiring `apply(self: *const P, z)` also keeps preconditioner application
structurally pure from the solver's perspective until a concrete mutable
preconditioner use case exists.

## 2026-04-05: Metric-aware Hodge star uses an explicit operator API while the existing operator context stays flat-specialized

**Decision:** Introduce metric parameterization at the Hodge-star layer via an
explicit `Metric(MeshType, mode)` type family plus `*_with_metric(...)` entry
points. The existing `hodge_star(...)`, `hodge_star_inverse(...)`, and
`OperatorContext` APIs remain the zero-cost flat specialization for now. The
Riemannian implementation covers interior-degree Whitney/Galerkin stars and the
top-degree diagonal star; primal 0-form Riemannian stars are deferred.

**Alternatives considered:**
1. Thread the metric through `OperatorContext`, codifferential, Laplacian, and
   all downstream examples immediately: rejected because it turns issue #85
   into a cross-component refactor spanning operator caching, assembled
   Laplacians, and example systems before the metric-aware star itself is
   validated.
2. Make the metric an optional runtime argument on the existing `hodge_star`
   surface: rejected because it obscures the flat fast path, weakens the type
   signal from `Metric(.flat)` vs `Metric(.riemannian)`, and makes the intended
   future generalization less explicit.
3. Force full Riemannian support for all degrees immediately: rejected because
   the current mesh stores only aggregate circumcentric dual volumes at
   vertices. A correct piecewise-metric `★₀` needs per-cell dual contributions
   or a redesigned dual-geometry representation, which is a separate topology
   issue.

**Rationale:** The vision says flatness must be an explicit specialization, not
an implicit assumption. The new metric-aware API makes that true at the Hodge
star boundary without paying a refactor tax across the rest of the operator
stack yet. Keeping the flat path unchanged preserves existing callers and
benchmarks. Deferring Riemannian `★₀` is honest: the missing data lives in the
mesh layer, so pretending otherwise would bake in a mathematically dubious
approximation.

## 2026-04-05: Boundary conditions are same-space projections; periodic uses explicit representatives

**Decision:** Boundary conditions in `src/operators/boundary_conditions.zig`
apply to a cochain and return the same cochain type
(`apply(allocator, input) -> InputType`). Dirichlet-like conditions derive
their constrained DOFs from topology-driven boundary masks. Periodic conditions
take an explicit representative map over DOFs instead of inferring pairings
from mesh geometry.

**Alternatives considered:**
1. In-place mutation only (`apply(*cochain)`): rejected because the BC library
   is easier to compose when each BC is a projection value rather than a
   mutation convention, and keeping a single public entrypoint avoids a
   parallel API.
2. Put automatic periodic pairing on `Mesh`: rejected because the topology
   layer owns incidence and geometry, not domain-identification policy.
   Baking periodic pairing into mesh construction would force rectangle/cuboid
   assumptions into a mesh-generic library.

**Rationale:** Dirichlet, PEC, and no-slip are all projections onto constrained
subspaces, so a same-space return type is the cleanest common model and fits
future BC composition. Periodic pairing is a simulation-setup choice rather
than a topological invariant. An explicit representative map keeps the
topology component honest and leaves room for later geometry-driven helpers
without freezing that policy into `Mesh` now.
