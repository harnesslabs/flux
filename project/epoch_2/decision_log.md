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
