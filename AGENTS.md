# flux

High-performance FEEC/DEC simulation framework in Zig. Targets accurate simulation of Maxwell's equations and incompressible Euler equations with exact discrete conservation laws.

See `project/initial.md` for the full architectural specification.

---

## Fast Orientation

Assume fresh context. Do **not** start by scanning the whole repository.

Read in this order:

1. This file for workflow rules and project expectations.
2. `project/components.md` for the codebase map and component boundaries.
3. Only the artifact that matches the current assignment:
   - implementation: the current issue/PR, then only the relevant component files
   - planning: the latest `project/epoch_*/roadmap.md`
   - design alignment: `project/vision.md` and `project/horizons.md`
   - decision logging: the current `project/epoch_N/decision_log.md`
4. Expand scope only when the current artifact proves you need more context.

Default rule: before any exploratory tool call, ask what the **smallest sufficient document or file set** is for the current task.

Do not read unrelated `src/` files "just to learn the repo". Learn the repo through `project/components.md`, then open only the component and direct dependencies that the task actually touches.

When designing APIs, prefer **one obvious way** to do something. Do not keep
parallel convenience and "real" APIs unless there is a concrete current need.
The project is pre-release; correctness and coherence matter more than
backward compatibility with transitional interfaces.

When extracting an abstraction from an issue or a duplicated example pattern,
do **not** stop at the nearest working seam. Pressure-test the candidate API
against the project's long-term target shape from `project/vision.md` and
`project/horizons.md` before committing to public nouns or ownership boundaries.
In practice:
- Prefer structural concepts over today's policy labels. A public type name
  should usually describe the underlying mathematical or computational role
  (`LinearSystem`, `EliminationMap`, `Integrator`) rather than one current
  use-case variant (`DirichletSystem`, `HeatRunner`, `PlaneScaffold`).
- Treat every new public type name as a mathematical claim. If it would feel
  wrong or embarrassingly narrow for a future example such as a mixed system,
  a high-dimensional phase-space problem, or a different solver family, the
  abstraction is probably still too shallow.
- Use example shrinkage as a quality signal, not the only success criterion.
  The stronger signal is that *new* examples would be straightforward to write
  without inventing new framework nouns.
- If an issue appears to call for a concrete local helper but the deeper
  structural abstraction is visible, prefer the deeper abstraction even if it
  requires a harder design pass up front. Local duplication removal is not
  enough by itself.
- Before finalizing an interface, ask explicitly: "What are the underlying
  concepts here? What would this look like for a PDE, discretization, or
  dimensionality we do not support yet?"

When an API is presented as degree-generic, dimension-generic, or otherwise
mathematically generic, pressure-test the implementation too — not just the
signature. A generic public API backed by a hardcoded case table over today's
supported dimensions/degrees is usually a smell. Prefer implementations derived
from the underlying structure or invariant (for example incidence, recursion,
or a comptime relation) unless there is a clear measured performance reason not
to.

---

## Hard Rules

- **Never merge a PR without explicit user approval.** After opening a PR,
  report the URL and stop. CI green is necessary but not sufficient — the user
  reviews every PR. Only run `gh pr merge` when the user explicitly says to
  (e.g., "merge it", "LGTM", "ship it"). This is non-negotiable.

---

## Build & Test

```sh
zig build                       # build
zig build test --summary all    # run all tests (with output)
zig build fmt                   # check formatting
zig build ci --summary all      # all CI checks: build + test + fmt
zig build run                   # run stub (directs to examples)
zig build run-maxwell           # run Maxwell examples (`--dim 2|3`)
zig build bench                 # run operator benchmarks (informational)
zig build bench -- --check      # compare against baselines, fail on >20% regression
zig build bench -- --update     # run benchmarks and save baselines.json
zig build bench -- --json       # output JSON to stdout
zig build bench -- --compare    # run same-run implementation comparisons
zig build docs                  # generate API docs to zig-out/docs/
zig build serve-docs            # build docs and serve at http://127.0.0.1:8080
```

### Benchmark Discipline

- Treat benchmark comparability as a first-class invariant. A perf claim is only valid if the compared benchmark rows measure the same thing with the same timing method.
- `bench/main.zig` uses two versioning layers:
  - `suite_version` is the JSON/baseline file format version.
  - Per-benchmark `version` is the measurement-method version for a single benchmark row.
- Bump a benchmark row's `version` when its measurement method changes: repetitions/batching, setup cost, workload size, warmup policy, data distribution, or benchmark semantics.
- Do **not** bump a benchmark row's `version` for an implementation-only change under the same measurement method. That is exactly when we want base-vs-PR comparisons to remain live.
- Never invalidate the whole suite because one row changed methodology. Preserve comparability for all unchanged rows.
- The default `zig build bench` suite should contain only surviving public library operations or surviving shipped example operations.
- Same-run implementation experiments belong behind `zig build bench -- --compare`, not in the default suite or baselines.
- When a perf PR adds a new benchmark with no base-branch counterpart, also include a same-run comparison that demonstrates the claimed win honestly, such as scalar-vs-default or old-path-vs-new-path within one run.
- PR benchmark comments must compare only rows with matching per-benchmark versions and must explicitly list skipped rows when versions differ.
- After a benchmark-method change merges, refresh `bench/baselines.json` on `main` with `zig build bench -- --update`.

---

## Vision Alignment

`project/vision.md` is the north star. Every design decision — discretization
choice, abstraction boundary, API shape — must be checked against it before
committing. If a decision deviates from the vision, it is not automatically
wrong, but it must be:

1. **Logged** in the epoch's `decision_log.md` as a deliberate, temporary move.
2. **Named**: state what vision principle it trades against and why the tradeoff
   is acceptable right now.
3. **Bounded**: note the condition under which it should be revisited.

A decision that conflicts with the vision and is not logged is a bug in the
project, not just the code.

`project/horizons.md` tracks validated but unscheduled architectural ideas.
Before finalizing any interface, check horizons to ensure the design does not
preclude a known future direction. When ideation produces a horizon-level
insight, add it there — not in the epoch roadmap.

Issue acceptance criteria are necessary but not sufficient design guidance.
They identify the concrete pain to remove, but they do not automatically define
the right public abstraction. If satisfying an issue literally would introduce
policy-shaped or example-shaped public nouns, step back and redesign around the
deeper reusable concept instead.

---

## Work Structure

Development has three levels of granularity. Each level has a definition of done.

### Sub-atomics: Commits

Commits are save states. Commit early, commit often. A commit should represent
a coherent, describable unit of progress — a passing test, a compiling refactor,
a new invariant — not a finished feature. They are local history: branchable,
cherry-pickable, revertable. Use them liberally.

Commit messages follow conventional commits (`feat:`, `fix:`, `test:`,
`refactor:`, `chore:`, `docs:`). The subject line says what changed; the body,
when needed, says why.

### Atoms: Issues

An issue is the smallest unit of releasable work. Each issue owns **3–5 tasks**
(checkboxes) that together deliver a complete capability. At minimum, one task
is a test task and one is an implementation task. Diagnostics, edge-case
handling, and documentation updates are bundled into the issue — not filed
separately.

**An issue is well-sized when:**
- It requires at least one non-obvious design choice (triggers `/decide`)
- It produces a PR that benefits from review (not trivially correct by inspection)
- It maps to roughly 5+ commits of meaningful progress

**An issue is too small when:**
- It's a single function addition with an obvious implementation
- The entire PR could be reviewed in under 2 minutes
- No design choices were made

**An issue is too big when:**
- It touches more than 2 component scopes (see `project/components.md`)
- The tasks list exceeds 7 items
- You can't state the single capability it delivers in one sentence

An issue is done when:

1. **The implementation exists** and compiles cleanly.
2. **Tests exist** that verify the invariant, not just the happy path. The
   testing hierarchy, in ascending order of proof strength:
   - Example-based unit tests (necessary but not sufficient)
   - Property-based / fuzz tests on random inputs (required for all operators)
   - `comptime` / type-level guarantees (promoted here whenever possible — a
     compile error is stronger than any runtime test)
3. **CI passes**: build, test, fmt, lint all green.
4. **The PR description is honest**: what was done, what was not done, what the
   known limitations are.

Open a draft PR when the branch is created, not when the work is done. This
keeps CI running continuously and creates a living dialogue thread for the work.

### Molecules: Pull Requests

Before moving a draft PR to ready-for-review, ask the downstream questions:

- Does any existing documentation need updating (including `project/vision.md`)?
- Does `README.md` need updating to stay consistent with the vision?
- Does `.github/settings.yml` (repo description, topics) need updating?
- Is there a new capability that warrants a new example demonstrating it
  end-to-end?
- Does this change the public API in a way that affects other modules?
- Does this introduce or resolve a decision that should be logged?
- Are there follow-on issues that should be opened now, while the context is
  fresh?

`README.md` and `.github/settings.yml` are public-facing artifacts that must
stay consistent with `project/vision.md`. When the vision evolves, they drift
unless explicitly checked.

A PR that introduces new atoms without answering these questions is incomplete,
even if CI is green. The molecule is not just the code — it is the code plus its
consequences in the broader system.

---

## Project Workflow

### Epochs → Milestones → Issues

- **Epoch**: A substantial body of work (~1 month or more). Living documents in `project/epoch_N/`. Tracked as a GitHub Project.
- **Milestone**: ~1–2 weeks of validated work. Each has an explicit mathematical acceptance criterion. 5–10 issues per milestone. Tracked as a GitHub Milestone.
- **Issue**: A complete capability with 3–5 tasks. Produces a PR with ~5+ commits. 5–10 per milestone. Tracked as GitHub Issues linked to the milestone.

A milestone is **done** when its acceptance criterion passes in CI — not when the code exists.

### Component Scopes

`project/components.md` maps each codebase component to its source files and
dependencies. Skills use this to limit context window usage. When working an
issue, load only the relevant component and its direct dependencies — do not
explore the full codebase.

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
  vision.md           # north star — design commitments and philosophy
  horizons.md         # validated but unscheduled architectural ideas
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
- In topology/operator code, prefer implementations phrased in the native discrete objects (`∂`, incidence, degree relations) over nested casework on `(dimension, degree)` when the former expresses the same rule directly

### Naming

- Full words, no abbreviations: `exterior_derivative`, not `ext_deriv`; `boundary_operator`, not `bop`
- Units and qualifiers appended in descending significance: `volume_dual_meters3`, `count_cells_max`

### Diagrams

- Use Mermaid (` ```mermaid ` code blocks) for all graphs, flowcharts, and dependency diagrams
- Never use ASCII box-drawing art — it breaks across fonts and renderers

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

### Development workflow

| Skill | Invoked by | When |
|-------|-----------|------|
| `/ideate` | user | Sharing raw ideas, brainstorming features, thinking aloud about project direction |
| `/epoch` | user | Planning upcoming work, discussing what to build next |
| `/milestone` | user | Epoch is planned, ready to create GitHub artifacts for a specific milestone |
| `/tackle` | user | Ready to start coding, asks "what should I work on" |
| `/decide` | agent (during `/tackle`) or user | A non-obvious architectural choice is made during implementation |
| `/review` | user | A PR is ready for review. Agent reviews implementation details; user reviews API + tests |
| `/retro` | user | Epoch boundary — reconstruct decisions, write retrospective, identify improvements |
| `/status` | user | Asks about progress, what's done, what's next, or what's blocking |

### Code quality (run on `main`, between epochs or mid-epoch)

| Skill | What it checks |
|-------|---------------|
| `/audit` | **Umbrella** — spawns all five lenses in parallel, consolidates results |
| `/audit-safety` | Assertion coverage, bounds, memory safety, numerical safety, invariant enforcement |
| `/audit-style` | Naming, dead code, stale references, comment quality, duplication |
| `/audit-cleanup` | API simplification, indirection removal, deduplication, shims, module/file cleanup |
| `/audit-perf` | Memory layout, allocation patterns, algorithmic efficiency, SIMD readiness |
| `/audit-tests` | Coverage gaps, property test quality, redundancy, convergence tests, comptime promotion |

### Workflow sequence

```
/ideate → /epoch → /milestone → /tackle (→ /decide as needed) → /review → /retro
                                                    ↑
                                          /audit-* (between epochs or mid-epoch)
```

**Proactive suggestions:**
- When the user shares loose ideas or "what if" thoughts → suggest `/ideate`
- After any planning discussion → suggest `/epoch` if no epoch doc exists yet
- After epoch doc exists and user asks about a milestone → suggest `/milestone`
- At the start of a coding session → suggest `/tackle` with a preview of the top issue
- After implementation work completes → suggest `/review` before merging
- At epoch boundaries → suggest `/retro` before starting next epoch planning
- When code feels messy or between epochs → suggest `/audit-*` to generate hygiene issues

---

## Target Directory Layout

```
src/
  root.zig          # library entry point — re-exports public API
  main.zig          # stub CLI (directs to examples)
  topology/         # CW complexes, incidence matrices
  forms/            # discrete k-forms / cochains
  operators/        # d, ★, and their compositions
  math/             # sparse linear algebra primitives
  io/               # VTK export
  integrators/      # generic time steppers (leapfrog, forward Euler)
  concepts/         # comptime concept validators
examples/
  maxwell/          # Maxwell examples (`root.zig` + README)
project/
  initial.md
  epoch_1/
  epoch_2/
  ...
```
