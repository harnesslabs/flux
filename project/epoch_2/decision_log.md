# Epoch 2 Decision Log

<!-- Append entries with /decide. Format: ## YYYY-MM-DD: <title> -->

## 2026-03-29: TimeStepper concept uses structural signature checking, not exact type equality

**Decision:** `TimeStepper(T)` validates `T.step` via `@typeInfo` — checking parameter
count, parameter types, and that the return is an error union with `void` payload — rather
than comparing `@TypeOf(T.step) == fn(Allocator, *State, f64) anyerror!void`.

**Alternatives considered:**
1. Exact type equality (`@TypeOf(T.step) == expected`): rejected because Zig distinguishes
   inferred error sets from `anyerror`. A function declared `fn(...) !void` gets a concrete
   error set, not `anyerror`, so no real implementation would match the exact type.
2. Duck typing (no concept check, just `anytype`): rejected because silent failures at
   call sites are harder to debug than a concept-level compile error.

**Rationale:** Structural checking accepts any error set in the return type, which is the
correct generality — an integrator's error set depends on which operators it calls. The
concept guarantees the parameter contract (allocator, state pointer, dt) while staying
permissive on the error set.
