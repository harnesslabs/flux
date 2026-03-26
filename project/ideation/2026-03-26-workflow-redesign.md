# Ideation: Workflow Redesign for Epoch 2

Date: 2026-03-26

## Context

Epoch 1 shipped all planned functionality (mesh, operators, Maxwell simulation) in under a week, but the development process revealed that issues were too small, most skills went unused, and the decision log stayed empty. This session re-examined every workflow skill and the issue/milestone granularity to prepare for epoch 2.

## Ideas explored

### Issue granularity: bigger issues, multiple tasks

**Layer:** constraint (process)
**Vision alignment:** aligned — "Tests as mathematical proof obligations" and the molecule checklist both benefit from issues substantial enough to generate meaningful review artifacts.
**Summary:** Epoch 1 issues were so small they were effectively single commits. Issues should own 3-5 tasks, bundle diagnostics with features, and require at least one non-obvious design choice. This makes `/tackle` → `/review` a meaningful pipeline instead of ceremony over trivial diffs.
**Key tension:** Sizing is partly feel-based, especially with AI writing code. Hard rules don't work — guidelines with invariants do.
**Verdict:** pursue
**Next action:** Update CLAUDE.md issue guidelines and `/epoch`/`/milestone` skill definitions with new sizing heuristics. An issue is well-sized when: (1) it requires at least one non-obvious design choice, (2) it produces a PR that benefits from review, (3) it maps to ~5+ commits. An issue is too big when it touches >2 component scopes, exceeds 7 tasks, or can't be stated in one sentence.

### `/tackle` redesign: test → stub → implement → review pipeline

**Layer:** constraint (process)
**Vision alignment:** aligned — directly enforces "write the mathematical property test before the implementation."
**Summary:** Formalize the flow: (1) open draft PR immediately, (2) write property tests encoding the acceptance criterion, (3) design public API as stubs, (4) implement until tests pass, (5) commit and push often for checkpointing. Draft PRs get CI running from the start. Frequent pushes create a safety net — bad turns can be reset without losing remote-tracked history. PRs should have ~5+ commits as a soft signal of adequate granularity.
**Key tension:** None — this is already close to what happened in epoch 1, it just wasn't formalized.
**Verdict:** pursue
**Next action:** Rewrite `/tackle` skill definition with this explicit flow.

### `/retro` skill for epoch retrospectives

**Layer:** constraint (process)
**Vision alignment:** aligned — "Every non-obvious architectural decision is logged" requires an enforcement mechanism. The retrospective is part of the epoch directory structure but had no skill to produce it.
**Summary:** `/retro` is user-invoked at epoch boundaries. It (1) reconstructs decisions from git history/PRs that were made but not logged, tagging them `[retroactive]`, (2) writes the epoch retrospective document, (3) identifies process improvements for the next epoch. Output feeds directly into `/epoch` for the next planning cycle.
**Key tension:** Retroactive decision logging is honest but second-class. The goal is to make `/decide` happen naturally during `/tackle` so `/retro` does less archaeology each time.
**Verdict:** pursue
**Next action:** Create `/retro` skill definition.

### Unused skills: `/milestone`, `/decide`, `/review`, `/status`

**Layer:** constraint (process)
**Vision alignment:** neutral
**Summary:** These skills weren't broken — they had nothing substantial to operate on. `/review` becomes valuable with meatier PRs. `/decide` should be invoked by the agent during `/tackle` when hitting design forks, not user-invoked. `/status` is useful when issues have multiple tasks and sessions span days. `/milestone` was unnecessary because milestones were too small to need separate planning.
**Key tension:** `/milestone` and `/status` may still not pull their weight. Keep them but don't force them.
**Verdict:** pursue (update, not remove)
**Next action:** Update `/decide` to be agent-invoked during `/tackle`. Keep `/status` and `/milestone` with updated granularity expectations. `/review` gets updated to match the user's review role (API shape + test coverage focus).

### `/audit` umbrella skill (safety, style, perf, tests)

**Layer:** constraint (process)
**Vision alignment:** aligned — "Safety > Performance > DX" maps to audit lenses.
**Summary:** Four separate skills under one namespace: `/audit safety`, `/audit style`, `/audit perf`, `/audit tests`. Each scans the codebase through its lens and produces findings that are either auto-fixable or issue-worthy. Run on `main` between epochs or mid-epoch. Separate from `/review` — review looks at a PR diff, audit looks at the whole codebase. Separate skills because each lens needs different context (saves context window).
**Key tension:** Four skills is a lot of surface area. But the alternative (one monolithic audit) wastes context loading irrelevant analysis.
**Verdict:** pursue
**Next action:** Create `/audit` skill definitions (4 sub-skills).

### Comptime trait system (#62)

**Layer:** architecture
**Vision alignment:** aligned — "comptime type safety for function spaces" and pluggable time integrators horizon.
**Summary:** Build a trait-like system using comptime interface checking and `@hasDecl`. This connects directly to the pluggable time integrators horizon already in `horizons.md`. The sequencing: build traits *during* Maxwell extraction to an example, not before — the concrete example validates the interface.
**Key tension:** None. This is the natural next step after epoch 1's concrete implementations.
**Verdict:** pursue
**Next action:** Candidate for epoch 2 planning. Sequence after or alongside Maxwell extraction.

### Code component scopes for agent context

**Layer:** UX (tooling)
**Vision alignment:** neutral (tooling concern)
**Summary:** Create `project/components.md` mapping each component to its source files, dependencies, and boundaries. Skills reference this to limit context window usage. When `/tackle` picks an issue in `domain/operators`, the agent loads only `src/operators/` + dependencies, not the whole codebase. The `domain/` GitHub labels already partially serve this purpose.
**Key tension:** This is a proxy for general context window management. The component map helps but isn't a complete solution.
**Verdict:** pursue
**Next action:** Create `project/components.md`. Update `/tackle` to read it and scope file access.

### Maxwell → example, not core module

**Layer:** architecture
**Vision alignment:** aligned — "Solver graph... users compose them. The framework does not prescribe the composition." Maxwell is a composition, not a framework primitive.
**Summary:** `src/em/` should eventually move to `examples/` as a demonstration of composing flux's operators. The framework provides the building blocks; Maxwell is an application. Leave in place for now, extract during epoch 2 when the trait system provides the right seams.
**Key tension:** None. This is clearly correct but depends on the trait system existing first.
**Verdict:** pursue (park for epoch 2 planning)
**Next action:** Epoch 2 candidate. Sequence after trait system.

### Compile-time Hodge decomposition

**Layer:** feature
**Vision alignment:** extends — pushes "comptime as a type system for mathematics" to decomposition-level structure.
**Summary:** The Hodge decomposition (ω = dα + δβ + γ) is a runtime computation. Comptime could tag results (e.g., `HarmonicForm` type) but can't compute the decomposition itself. Comptime Betti numbers would require comptime meshes — a large architectural commitment with narrow utility. The decomposition may already be implicit in the FEEC framework. Parked as a "stoner thought" — not concrete enough to act on.
**Key tension:** No clear problem this solves today.
**Verdict:** park
**Next action:** None. Revisit if a concrete simulation needs to branch on harmonic subspace structure.

### Stale issue references in codebase

**Layer:** constraint (hygiene)
**Vision alignment:** neutral
**Summary:** Inline comments reference issue numbers and possibly other artifacts that may be stale. This is a cleanup task, not a design decision.
**Verdict:** park
**Next action:** `/audit style` would catch this. No separate action needed.

## Connections

- **Granularity drives everything.** Bigger issues → meaningful `/decide` during `/tackle` → useful `/review` → honest `/retro`. The skill pipeline works when the atoms are the right size.
- **Trait system + Maxwell extraction** are two sides of the same coin. The traits define the interface; Maxwell extraction validates it. Sequence them together in epoch 2.
- **Component scopes + `/audit` sub-skills** both address context window management from different angles. Components limit what the agent reads; audit sub-skills limit what analysis is loaded.
- **`/tackle` pipeline + `/review` role clarity** form a contract between agent and user: the agent owns implementation quality, the user owns API shape and test design review.

## Open questions

- What is the right epoch duration? Epoch 1 was planned for ~1 month and completed in <1 week. Should epoch 2 be more ambitious, or should we accept that AI-assisted velocity means shorter epochs?
- Should `/milestone` remain a separate skill, or should milestone creation fold into `/epoch`?
- How do we measure whether the new issue sizing guidelines are working? Is "the decision log has entries" a sufficient signal?
