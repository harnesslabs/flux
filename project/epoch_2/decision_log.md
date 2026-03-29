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

## 2026-03-29: OperatorConcept validates types, not functions

**Decision:** `OperatorConcept(Op)` validates a *type* (namespace struct) with `InputType`,
`OutputType`, and `apply(Allocator, InputType) !OutputType`. Existing free-function operators
are wrapped into concept-conforming types via parameterized constructors (e.g.,
`ExteriorDerivativeOperator(MeshType, k, Duality)`). The free functions remain the
implementation layer; wrappers are the generic interface layer.

**Alternatives considered:**
1. Validate function signatures directly (concept over functions, not types): rejected
   because functions lack metadata. A function `fn(Allocator, anytype) !U` doesn't carry
   its input/output types as inspectable declarations. The type approach gives generic code
   access to `Op.InputType` and `Op.OutputType` for composability checking.
2. Refactor existing operators into types directly: rejected as too invasive for this PR.
   The existing `anytype`-based free functions work well with `compose.chain()`. The
   wrapper approach adds the type-level interface without disrupting existing code.
3. No wrapper, concept only: would leave the concept without any conforming types in the
   codebase, making it untested against real operators.

**Rationale:** Types carry metadata that functions don't. Generic infrastructure (solver
graphs, operator algebras) needs to inspect input/output types at comptime. The wrapper
pattern mirrors `MaxwellLeapfrog(MeshType)` wrapping `leapfrog_step` for `TimeStepStrategy`.

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

**Decision:** Comptime concept validators use the `*Concept` suffix: `OperatorConcept`,
`MeshConcept`. This is slightly inconsistent with the existing `TimeStepStrategy` name
(which serves the same role) but more descriptive.

**Alternatives considered:**
1. `*Strategy` suffix for all: rejected because "strategy" implies a behavioral pattern
   (Strategy Pattern), not a type contract. "Concept" is the standard C++ term for
   compile-time type constraints and maps well to Zig's comptime validation.
2. Rename `TimeStepStrategy` → `TimeStepConcept` now: deferred to avoid scope creep.
   Should be done in a follow-up issue to unify naming.

**Rationale:** `*Concept` is the clearest name for "comptime function that validates a
type satisfies an interface." Unifying `TimeStepStrategy` is tracked as a follow-up.
