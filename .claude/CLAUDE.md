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

## GitHub Labels

Labels are defined in `.github/settings.yml` (managed by Probot) and auto-applied to PRs via `.github/labeler.yml` based on changed file paths.

| Prefix | Purpose | Values |
|--------|---------|--------|
| `type/` | Kind of work | `bug`, `feature`, `invariant`, `perf`, `docs`, `refactor`, `test`, `ci` |
| `phase/` | Roadmap phase | `1` through `5` |
| `domain/` | Codebase area | `topology`, `forms`, `operators`, `io`, `em`, `fluid`, `build` |
| `priority/` | Urgency | `high`, `medium`, `low` |
| `status/` | Workflow state | `blocked`, `needs-decision`, `good-first-issue` |

When creating issues, always set `type/`, `phase/`, `domain/`, and `priority/`. Use labels for filtering:
```sh
gh issue list --label "phase/1,type/feature"   # current phase features
gh issue list --label "type/invariant"          # all math property work
gh issue list --label "status/blocked"          # what's stuck
```

## Issue Templates

Three templates in `.github/ISSUE_TEMPLATE/`:
- **implementation** — new functionality; requires a concrete acceptance criterion
- **bug** — broken/incorrect behavior; asks for the violated invariant
- **invariant** — mathematical property to implement or verify; requires a test plan

Blank issues are disabled. Every issue must use a template.

## CI

Four jobs run in parallel on every push and PR (`.github/workflows/ci.yml`):

| Job | Command | What it catches |
|-----|---------|-----------------|
| `build` | `zig build` | Type errors, comptime failures |
| `test` | `zig build test` | All `test` blocks |
| `fmt` | `zig fmt --check src/ build.zig` | Formatting |
| `lint` | `zig ast-check` on all `.zig` files | Parse/AST errors without full codegen |

All four must pass. Branch protection is configured in `.github/settings.yml`.

---

## Skills

Six slash commands drive the development workflow. Suggest the appropriate one proactively when the context calls for it — don't wait to be asked.

| Skill | When to suggest |
|-------|----------------|
| `/epoch` | User wants to plan upcoming work, discusses what to build next, or asks about direction |
| `/milestone` | Epoch is planned and user is ready to create GitHub artifacts for a specific milestone |
| `/tackle` | User wants to start coding, asks "what should I work on", or is ready to close an issue |
| `/review` | A PR is open, user asks about code quality, or a branch has been pushed |
| `/status` | User asks about progress, what's done, what's next, or what's blocking |
| `/decide` | A non-obvious architectural choice is made during any session |

**Proactive suggestions:**
- After any planning discussion → suggest `/epoch` if no epoch doc exists yet
- After epoch doc exists and user asks about a milestone → suggest `/milestone`
- At the start of a coding session → suggest `/tackle` with a preview of the top issue
- After implementation work completes → suggest `/review` before merging
- When an architectural choice comes up mid-session → suggest `/decide` immediately, don't let it slide

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
