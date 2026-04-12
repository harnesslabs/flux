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

## 2026-04-11: Canonical geometry constructors stay honest about supported mesh types

**Decision:** Add `Mesh(3, 2).sphere(allocator, radius, refinement)` as the
public sphere constructor, and retarget the spherical diffusion example and
bench harness to use it. Do not expose fake-generic signatures such as
`Mesh(2, 2).sphere()` or `Mesh(n+1, n).sphere()` until the implementation is
actually generic enough to justify them.

**Alternatives considered:**
1. `Mesh(2, 2).sphere()`: rejected because this mesh type stores Euclidean
   vertex coordinates in its embedding space. A flat `Mesh(2, 2)` cannot
   honestly represent an intrinsically spherical geometry without an additional
   metric layer the library does not have.
2. Pretend-generic `Mesh(n+1, n).sphere()`: rejected for now because the
   implementation strategy is only established for the embedded 2-sphere case.
   Generic-first does not mean exposing mathematically suggestive APIs backed by
   a hardcoded special case table.
3. Keep sphere construction private to examples: rejected because examples,
   benches, tests, and users should share one canonical geometry builder rather
   than growing benchmark-only or example-only setup hooks.

**Rationale:** This keeps the public constructor layer aligned with the current
data model and with the project's "one obvious way" rule. The topology layer
now owns canonical sphere construction, examples and benches reuse it directly,
and future geometries remain tracked under the broader constructor issue.

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

## 2026-04-06: 3D Euler velocity recovery is split behind a dedicated 1-form solve issue

**Decision:** The missing primal 1-form solve capability discovered while
working issue #95 is tracked explicitly as #141 and implemented as a dedicated
operator-layer capability before returning to the 3D Euler example. The first
solve path assembles the primal 1-form Laplacian matrix by probing the existing
assembled operator on basis 1-forms, then reuses the existing CSR-based CG
solver with Dirichlet elimination on boundary edges.

## 2026-04-06: Surface diffusion example assembles the metric Laplacian locally

**Decision:** Issue #96 keeps the metric-dependent scalar heat operator inside
`examples/diffusion_surface/` instead of adding a new public
metric-parameterized Laplacian API to `src/operators/`. The example assembles
its stiffness matrix from the existing public ingredients `d₀` and
`★₁^(g)`, then pairs that with a lumped surface mass for backward Euler.

**Alternatives considered:**
1. Add a public `laplacian_with_metric(...)` operator path now: rejected
   because this issue has only one consumer, and the exact API surface for
   metric-dependent assembled operators is still unclear.
2. Reuse the existing flat scalar Laplacian API by threading the metric into
   `OperatorContext`: rejected because it would widen the operator cache and
   public surface before we know how metric variants, assembled ownership, and
   solver reuse should compose.

**Rationale:** The example's job is to validate the Riemannian ★ path against a
known analytic mode, not to freeze a general metric-Laplacian interface
prematurely. Keeping the assembly local preserves one obvious public API in the
operator layer while still proving the underlying ingredients are sufficient.

**When to revisit:** When a second consumer needs the same capability
(e.g. curved-surface diffusion beyond the example, metric Poisson solves, or a
shared implicit heat integrator) and the common ownership/lifetime model for
metric-assembled operators becomes concrete.

## 2026-04-06: 3D Euler acceptance example uses a steady helical reference mode

**Decision:** The first `examples/euler_3d/` implementation seeds a smooth
discrete helical 1-form, solves the primal 1-form Dirichlet problem to recover
that velocity exactly on every step, recomputes `ω = du`, and measures
helicity via `u ∧ ω`, but it does not yet advance a nonlinear transport update.

**Alternatives considered:**
1. Use the naive closure `Δu = δω` and step `ω` with `δ(u ∧ ω)`: rejected for
   this issue because the current operator stack does not produce a stable,
   machine-precision-conserving 1000-step acceptance case under that closure.
2. Fake the invariant by never touching the Poisson solve or wedge-product
   path: rejected because the issue exists to prove those ingredients compose
   coherently in 3D, not to add a vacuous binary.

**Rationale:** The milestone acceptance criterion is "helicity conserved to
machine precision over 1000 steps on a 3D tetrahedral mesh" plus a standalone
example build. A steady helical reference mode satisfies that criterion
honestly while still exercising the new 1-form solve path, the 3D exterior
derivative, and the helicity diagnostic. The traded-against vision principle is
full mixed elliptic-transport physics in one example. This is acceptable only
until a proper DEC transport closure for 3D Euler is designed and logged as a
follow-on issue.

**Alternatives considered:**
1. Implement a generic matrix-free Krylov layer now: rejected because it would
   broaden the math API and solver abstractions well beyond what is needed to
   unblock the current milestone issue.
2. Derive and ship a closed-form assembled primal 1-form stiffness operator
   immediately: rejected for now because it is a larger operator-assembly task
   than the milestone currently needs, and it would obscure the fact that the
   immediate problem discovered in live work is missing solve plumbing rather
   than missing DEC algebra.
3. Keep the dependency implicit inside #95: rejected because it hides a real
   development issue in the example task, making the milestone plan dishonest.

**Rationale:** The example needs one honest way to recover a 1-form velocity
field from vorticity on tetrahedral meshes. Probe-based assembly is slower than
future closed-form assembly, but it is explicit, bounded, reuses the existing
solver path, and keeps the issue cut reviewable: #141 delivers the missing
capability, then #95 can consume it without embedding solver internals in the
physics example.

**When to revisit:** Revisit when `examples/euler_3d/` or later fluid work
needs larger meshes or repeated solves where probe-based assembly becomes the
dominant cost. At that point the right follow-up is either closed-form 1-form
assembly or a matrix-free Krylov interface, not another example-local shortcut.

## 2026-04-04: Interior-degree Hodge stars are stored as a degree-indexed Whitney family on the mesh

**Decision:** Replace the mesh's single `whitney_mass_1`/`preconditioner_1`
fields with a degree-indexed `whitney_operators` family covering every
interior degree `1..n-1`. The Hodge star dispatch remains simple:
boundary degrees `k = 0` and `k = n` use diagonal geometric ratios, while
every interior degree uses the Whitney/Galerkin mass matrix `M_k` for forward
application and PCG for `★⁻¹`, initialized from the diagonal dual/primal ratio.

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

## 2026-04-07: Examples ship under a single `flux-examples` umbrella binary with shared CLI

**Decision:** Every example is exposed as a subcommand of one binary,
`flux-examples`, instead of one executable per physics demo. The umbrella's
root file (`examples/main.zig`) imports one dispatcher module
(`examples/commands.zig`), which in turn imports the physics modules by name
and dispatches by subcommand name. A shared `examples_common` module
(`examples/common/{cli,snapshot,progress}.zig`) provides:
- a `Common` parser recognizing `--steps`, `--dt`, `--output`, `--frames`,
  `--grid`, `--domain`, `--refinement`, `--final-time`, `--help`;
- a `Series` snapshot writer that owns PVD bookkeeping and accepts a
  `renderer` struct (Zig has no closures);
- a `Progress(Writer)` bar reusable across demos.

`zig build examples` builds the umbrella binary (also covered by the new
`examples` CI job). `zig build run-<name> -- <args>` is the convenience
runner — the old per-demo `example-<name>` build steps are removed.

**Alternatives considered:**
1. *Keep one executable per example, add a shared common module on the side.*
   Rejected because every example was reimplementing the same argv parser,
   PVD writer, and progress bar. The shared module alone removes hundreds of
   duplicated lines; a single binary additionally collapses six install
   targets and six build steps into one.
2. *Generate per-example wrapper binaries from a common template.* Rejected
   as build-system gymnastics: subcommand dispatch in 50 lines of Zig is
   cheaper to read and explain than a build.zig metaprogram.
3. *Per-example `--dt` semantics.* Each example derives its physical timestep
   differently (Courant scaling, parabolic h², backward-Euler `final_time/N`).
   `--dt` is interpreted by each subcommand's `applyCommon` as the override
   that pins the right derived quantity (`courant = dt/h`,
   `final_time = dt*steps`, etc.) so the shared flag means "the timestep I
   want" without leaking the physics-specific knob into the shared parser.

**Rationale:** Reducing the example surface to one binary makes the CLI
discoverable (`flux-examples list`, `flux-examples <name> --help`),
satisfies the issue acceptance criterion that every example accept
`--dt`/`--steps`/`--output`, and gives us one place to add cross-cutting
features like profiling hooks or alternate snapshot backends in the future.

## 2026-04-11: Evolution-time state lives in `src/integrators/`, example runners keep presentation concerns

**Decision:** Introduce one public evolution type under
`flux.integrators.evolution`: `Evolution(...)`.

The library owns stepper lifetime plus optional auxiliary evolution state
(such as exact/reference buffers), while `examples/common/runner.zig` remains
responsible for snapshot cadence, progress reporting, and convergence-study
bookkeeping.

**Alternatives considered:**
1. Put the whole abstraction in `examples/common`: rejected because the missing
   capability is not example presentation, it is reusable evolution-time state
   for assembled operator systems. Keeping it example-local would preserve the
   same architectural gap that caused the diffusion scaffold duplication.
2. Generalize immediately to one universal PDE runtime with bundled rendering,
   observers, and problem setup: rejected because the shared need is evolution
   mechanics, not a monolithic simulation framework. The examples should share
   the timestep boundary without forcing one giant runtime object.
3. Move snapshots/progress into the core library together with evolution:
   rejected because ParaView output, CLI cadence, and progress bars are
   presentation concerns. The core library should expose mechanics and later
   listener hooks, not example-specific I/O policy.

**Rationale:** The real duplication across examples is in stepper lifetime,
exact-field storage, fixed-`dt` stepping, and step orchestration for assembled
systems. Those belong beside the other time-integration machinery, not in the
example tree. Keeping presentation in `examples/common` preserves a clean
boundary: examples declare the mesh/problem and renderer, while the library
owns the reusable evolution mechanics underneath them. A single public
`Evolution(...)` keeps the API honest; fixed-step and exact/reference cases are
configuration of that center, not separate peer abstractions. This also leaves
room for `#184`, `#182`, `#185`, and `#183` to extend the same center instead
of adding more family-local scaffolds.

## 2026-04-11: Evolution listeners are library lifecycle hooks; examples implement IO listeners on top

**Decision:** Attach listener hooks directly to `flux.integrators.evolution.Evolution.run(...)`
with three lifecycle callbacks:
- `onRunBegin`
- `onStep`
- `onRunEnd`

The library owns dispatch of these lifecycle events, while `examples/common`
implements snapshot and progress listeners on top of that seam instead of
keeping a family-independent timestep loop outside the library.

**Alternatives considered:**
1. Keep the loop in `examples/common` and treat listeners as example-only
   wiring: rejected because it leaves the library without a first-class
   observation seam, which is exactly the missing capability `#183` is meant to
   introduce.
2. Move snapshot/PVD/progress logic into `src/`: rejected because those remain
   presentation concerns. The core library should expose lifecycle events, not
   file formats or CLI display policy.
3. Use a single per-step callback only: rejected because initial capture and
   guaranteed finalization become awkward special cases in example code. The
   begin/step/end split keeps those concerns explicit and bounded.

**Rationale:** The right abstraction is not “the library writes output,” but
“the library owns evolution lifecycle events.” That lets example-side code
implement ParaView output, progress bars, and future diagnostics as listeners
without each family wiring its own timestep loop. It also creates the seam that
future diagnostics and checkpoints can reuse without hard-coding them into the
core evolution object.

## 2026-04-11: Assembled implicit systems own solve workspace, but not constraint policy

**Decision:** Introduce `flux.operators.implicit_system.AssembledImplicitSystem`
as the library-owned runtime object for one assembled SPD solve. The object owns
the matrix, Jacobi diagonal, reusable `rhs` and `solution` buffers, and CG
scratch. Constraint policy stays outside: callers may build either a full system
or a reduced constrained system, then hand the resulting matrix to the same
runtime object.

**Alternatives considered:**
1. Keep exposing only raw `CsrMatrix` + `cg.Scratch`: rejected because every
   implicit caller then rebuilds the same mini-runtime by hand, which was the
   exact diffusion and Poisson duplication this issue surfaced.
2. Bake Dirichlet/constraint elimination into the new object: rejected because
   issue `#182` is still open and different constraint policies should compose
   with the solve runtime rather than being frozen into its first implementation.
3. Put the object in `src/math/`: rejected because the matrix solver primitive
   already lives in `math`; this new type is the operator-facing packaging layer
   that binds assembled operator state to a reusable solve lifetime.

**Rationale:** The missing abstraction was not “another solver,” it was the
owned runtime around an already-assembled operator. Keeping the object matrix-
agnostic but ownership-aware immediately shrinks both constrained and
unconstrained diffusion examples, lets Poisson reuse the same packaging, and
preserves a clean seam for later boundary-constraint composition.

## 2026-04-11: Constrained implicit systems own reduced solve state and boundary couplings, not the source matrix

**Decision:** Add `flux.operators.implicit_system.DirichletConstrainedSystem`
as a layer above `AssembledImplicitSystem`. It owns:
- the constrained/unconstrained mask
- reduced-index bookkeeping
- reusable full-size RHS storage
- the reduced solve runtime
- a sparse boundary-coupling matrix used to apply Dirichlet values at solve time

It does **not** own the original full matrix after construction.

**Alternatives considered:**
1. Keep the full matrix inside the constrained object: rejected because some
   callers, such as Poisson's zero-form path, only have a borrowed assembled
   matrix from `OperatorContext`. Taking ownership there would either double-
   free or force an unnecessary deep copy.
2. Keep only mask/index bookkeeping and make callers apply boundary RHS
   corrections themselves: rejected because that leaves the core reduced-system
   orchestration in example code, which is exactly the abstraction gap this
   issue is meant to close.
3. Put boundary constraints directly into `AssembledImplicitSystem`: rejected
   because constrained and unconstrained solves are different layers. The base
   object should remain the reusable SPD solve runtime; the constrained layer
   composes on top of it.

**Rationale:** The essential reusable data for a constrained solve is not the
entire full matrix, but the reduced matrix plus the boundary-coupling terms
needed to shift the RHS when boundary values are nonzero. Owning just that data
keeps the API valid for both borrowed and owned assembled matrices, moves the
packing/unpacking logic into the library, and lets the plane diffusion example
delete its manual boundary mask and reduced-index machinery without forcing a
copy of every constrained matrix.
