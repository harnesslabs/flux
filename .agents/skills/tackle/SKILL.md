---
name: tackle
description: Select the highest-priority unblocked issue for the current milestone, work it to a reviewable PR, and stop for review.
---

Find and work the highest-priority open issue to completion — one issue, one branch, one PR, CI passes, stop for review.

## Find the issue

1. Read `project/components.md` to understand codebase scope boundaries.
2. Read `project/epoch_*/roadmap.md` for the current epoch to understand milestone priorities and acceptance criteria.
3. Read `project/horizons.md` — before writing any interface, check that it does not preclude a known future direction.
4. Read `project/vision.md` when the work introduces or reshapes a public abstraction. Treat issue wording as the local symptom, not automatically as the correct abstraction boundary.
5. Get live state:
   ```sh
   gh milestone list --state open
   gh issue list --state open --label "priority/high" --json number,title,labels,milestone
   gh issue list --state open --label "status/blocked" --json number,title,labels,milestone
   ```
6. Select the best issue to work next using this priority order:
   - On the current milestone (check which milestone is earliest in the roadmap)
   - `type/invariant` preferred — these are on the critical path for acceptance criteria
   - `priority/high` over `priority/medium`
   - Issues that unblock other `status/blocked` issues
   - Skip `status/blocked` issues (they can't be worked until unblocked)
   - Non-epoch bugs (`type/bug`) may take priority if they are breaking something

7. Present the chosen issue to the user: number, title, acceptance criterion, and a one-sentence rationale for why this one. Ask for confirmation or redirection.

## Scope the work

On confirmation:

1. Read the full issue body: `gh issue view <number>`
2. Identify the **component scope** from `project/components.md` based on the issue's `domain/` label. Name the component and direct dependencies explicitly before opening source files.
3. State the deeper structural concept you expect the issue to touch. If the issue title/body names a policy or one example family, ask whether that is the real abstraction or only the current manifestation.
4. Create a branch and immediately open a draft PR:
   ```sh
   git checkout main && git pull
   git checkout -b <number>-<short-slug>
   git commit --allow-empty -m "chore: open draft PR for #<number>"
   git push -u origin <number>-<short-slug>
   gh pr create --draft \
     --title "<imperative title>" \
     --body "Closes #<number>

   ## What

   <one sentence: what this PR delivers>

   ## Acceptance criterion

   <from the issue>

   ## Tasks

   - [ ] Write property tests encoding the acceptance criterion
   - [ ] Design public API (stubs)
   - [ ] Implement
   - [ ] CI green
   " \
     --label "<type/X,domain/X,priority/X>"
   ```
   Slug format: lowercase, hyphens, ≤5 words. Example: `42-csr-incidence-matrix`.

5. Read all relevant source files within the component scope before writing any code.
   Start with the files named in the issue and the component entry in `project/components.md`. Only expand to additional direct dependencies when the current files prove they are needed. Do not read unrelated modules "for context".

## Phase 1: Tests

Write the tests first. These define what "done" means.

1. Identify the mathematical invariant or behavioral property from the acceptance criterion.
2. Write property tests that verify this invariant. Use random inputs where applicable.
3. Write example-based tests for edge cases (empty mesh, boundary cells, degenerate inputs).
4. Tests should fail at this point — the implementation doesn't exist yet.
5. Commit and push:
   ```sh
   git add <test files>
   git commit -m "test(domain): add failing tests for <invariant>"
   git push
   ```

## Phase 2: API design

Design the public interface before implementing internals.

1. Write function signatures, struct definitions, and type-level constraints as stubs.
2. Use `@compileError("not yet implemented")` or `unreachable` for function bodies.
3. Verify the tests compile against the stubs (they should fail at runtime, not compile time).
4. Pressure-test every new public noun against at least one uncomfortable future case from `project/vision.md` / `project/horizons.md`:
   - a different PDE family
   - a different dimensionality
   - a different constraint/policy choice
   - a solver-graph composition use case
   If the name or ownership model would obviously need replacement there, redesign before implementing.
5. **Default-forward checkpoint:** present the API surface to the user briefly, but continue without waiting for confirmation unless the design is materially uncertain, introduces a cross-component interface, or would be expensive to unwind if wrong. Ask for confirmation only in those cases.
6. If a non-obvious design choice was made (struct layout, ownership model, comptime parameter choice), log it immediately:
   - Read the current epoch's `decision_log.md`
   - Append the decision in the standard format (see `/decide`)
   - Commit the decision log alongside the stubs
7. Commit and push:
   ```sh
   git add <source files> <decision log if updated>
   git commit -m "feat(domain): stub API for <capability>"
   git push
   ```

## Phase 3: Implement

Fill in the stubs until tests pass.

1. Implement in small steps. Each step must compile and pass all existing tests before moving on.
2. **Commit after every coherent unit of progress** — a working constructor, a helper function, a passing test case. Target roughly 50–150 lines per commit; if you are approaching 200+ lines uncommitted, stop and commit what you have.
3. `zig build ci --summary all` between every commit. No exceptions.
4. Push after every 1–2 commits to maintain remote checkpoints.
5. If you realize an approach is wrong, you can reset to an earlier commit — the remote has the history. State what you're resetting and why before doing it.
6. When all acceptance criterion tests pass, that is the final implementation commit.
7. If additional non-obvious decisions were made during implementation, log them to the decision log.
8. If the issue is performance work, benchmark discipline is part of the implementation:
   - Keep at least one honest comparable benchmark row for every claimed speedup when possible.
   - If benchmark methodology changes for a row, bump that row's per-benchmark `version` only for the affected benchmarks.
   - Do not claim a base-vs-PR speedup from a row whose benchmark `version` changed.
   - If the base branch lacks a benchmark for the claimed win, add a same-run comparison benchmark so the PR still demonstrates the improvement.
9. Maintain a short list of **scoped-out follow-on work** discovered while implementing:
   - API annoyances, missing topology hooks, duplicated logic, naming drift, or cleanup that became obvious only because you were in the code
   - Do **not** expand the current PR to include these unless they are required to finish the issue correctly
   - Before creating anything, check whether the follow-on is already tracked
   - If it is real, untracked, and worth doing later, open a **non-milestone** GitHub issue with the normal labels (`type/`, `domain/`, `priority/`) and explicitly leave milestone unset

10. If the implementation solves the issue by introducing a narrow policy-shaped abstraction, stop and redesign before finalizing. The goal is to remove the symptom by exposing the deeper reusable concept, not by wrapping today's special case in a new public type.

The rhythm is: **write → test → commit → push → repeat**.

## Phase 4: Finalize

1. Run `zig build ci --summary all` one final time.
2. Update the draft PR description:
   - Check all completed task boxes
   - Add a "## Decisions" section if any were logged
   - Add a "## Limitations" section for known gaps
3. Answer the molecule checklist:
   - Does this introduce an interface that conflicts with a known horizon in `project/horizons.md`?
   - Does existing documentation need updating?
   - Does this change the public API in a way that affects other modules?
   - Are there follow-on issues that should be opened?
4. If you discovered real follow-on work during implementation and it is not already tracked, open the non-milestone issue(s) now and mention them in the PR description.
5. Mark the PR as ready for review:
   ```sh
   gh pr ready <number>
   ```
6. Report the PR URL and **STOP**. Do not merge.

The user reviews API shape and test coverage. The agent (via `/review`) reviews implementation details.

## After review

Merging is the user's decision. Only merge when the user explicitly says to (e.g., "merge it", "LGTM", "ship it"). Then:

1. Check CI: `gh pr checks <number>` — all checks must pass.
2. Merge:
    ```sh
    gh pr merge <number> --squash --delete-branch
    ```
3. Confirm the issue was auto-closed: `gh issue view <number>`
4. Report: what was merged, what the next recommended issue is.

## Constraints

- **Never merge without explicit user approval.** CI green is necessary but not sufficient.
- Never skip `zig build ci --summary all` before pushing.
- Never merge with failing CI.
- Never put more than one issue's work in a single PR.
- Never write more than ~200 lines without compiling, testing, and committing.
- Stay within the component scope identified at the start. If you find yourself needing to read or modify files outside the scope, state why and get confirmation.
- Draft PR opens immediately. Commits push frequently. The remote is your checkpoint system.
- A PR with fewer than 3 commits is a signal that either the work was trivial (issue was too small) or the agent was not committing incrementally. Neither is acceptable for well-sized issues.
