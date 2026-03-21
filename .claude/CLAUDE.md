# flux

High-performance FEEC/DEC simulation framework in Zig. Targets accurate simulation of Maxwell's equations and incompressible Euler equations with exact discrete conservation laws.

See `project/initial.md` for the full architectural specification.

---

## Build & Test

```sh
zig build          # build
zig build test     # run all tests
zig build run      # run executable
```

---

## Project Workflow

### Epochs → Milestones → Issues

- **Epoch**: ~1 month of work. Living documents in `project/epoch_N/`. Tracked as a GitHub Project.
- **Milestone**: ~1 week of validated work. Each has an explicit mathematical acceptance criterion. Tracked as a GitHub Milestone.
- **Issue**: 2–6 hours of focused work. 10–20 per milestone. Tracked as GitHub Issues linked to the milestone.

A milestone is **done** when its acceptance criterion passes in CI — not when the code exists.

### Decision Log

Every non-obvious architectural decision goes in `project/epoch_N/decision_log.md` with:
- The decision made
- Alternatives considered
- Rationale at the time

This is a first-class artifact. Before refactoring something that seems wrong, check the decision log — there may be a reason.

### Epoch Directory Structure

```
project/
  initial.md          # original architectural specification (read-only reference)
  epoch_1/
    roadmap.md        # milestones, goals, acceptance criteria
    decision_log.md   # architectural decisions + rationale
    retrospective.md  # written at epoch end — what held, what didn't
  epoch_2/
    ...
```

---

## Mathematical Invariants

These must hold exactly (to machine precision where applicable) and must be verified by tests:

- **dd = 0**: applying the exterior derivative twice yields zero (cohomological identity); test on random inputs
- **∇·B = 0**: discrete magnetic divergence identically zero — enforced structurally, not approximately
- **Circulation conservation**: total circulation over macroscopic loops preserved for incompressible Euler flows

Property-based tests on random meshes and random cochain inputs are the primary correctness mechanism for discrete operators. An operator is not implemented until this kind of test exists and passes.

---

## Code Conventions

### Zig

- Use `comptime` to enforce k-form type safety — a function accepting a 1-form must reject a 2-form at compile time, not runtime
- SoA layout via `std.MultiArrayList` for all mesh entities; justify explicitly if deviating
- Explicit allocators everywhere — no hidden heap allocations
- `const` by default; mutability is the exception and should be obvious
- Assert preconditions and invariants with `std.debug.assert`; assertions are not optional
- All loops must be bounded — no unbounded iteration

### Naming

- Full words, no abbreviations: `exterior_derivative`, not `ext_deriv`; `boundary_operator`, not `bop`
- Units and qualifiers appended in descending significance: `volume_dual_meters3`, `count_cells_max`

### Comments

- Explain *why*, not what — the code says what
- Mathematical content uses Unicode symbols: ∂, ∇, ∈, ★, ∧, ⟨⟩, Ω
- Cite sources with author, title, and equation number — do not cite fabricated URLs

### Testing

- Write the mathematical property test before the implementation
- Use `std.testing.expectEqual` for exact checks; tolerance-based comparison for floating-point with explicit epsilon
- Test names should state the invariant: `test "dd = 0 for random 1-forms on triangular mesh"`

---

## Target Directory Layout

```
src/
  root.zig          # library entry point — re-exports public API
  main.zig          # CLI entry point
  topology/         # CW complexes, incidence matrices (Phase 1.1)
  forms/            # discrete k-forms / cochains (Phase 1.2)
  operators/        # d, ★, and their compositions (Phase 1.3–1.4)
  io/               # VTK export (Phase 2.2)
  em/               # electromagnetics integrators (Phase 3)
  fluid/            # fluid dynamics integrators (Phase 4)
project/
  initial.md
  epoch_1/
  epoch_2/
  ...
```
