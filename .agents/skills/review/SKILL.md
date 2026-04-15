---
name: review
description: Review a pull request with emphasis on implementation details, mathematical correctness, type safety, invariants, and scope discipline.
---

Review a PR for the flux project with mathematical and code correctness scrutiny.

Usage: /review <PR number>

## Context: division of review labor

The **user** reviews API shape and test coverage — they check that the interface is clean and the tests are comprehensive. They rarely need to read implementation blocks in detail.

The **agent** (this skill) reviews implementation details with a fine-tooth comb — catching subtleties the user will miss at a glance. This is the agent's primary value in review.

Structure your review accordingly: lead with the implementation-detail findings that the user is unlikely to catch on their own. Don't spend the review restating things visible from the PR description.

## Gather context

1. Fetch PR metadata and changed file paths first:
   ```sh
   gh pr view <number> --json number,title,body,labels,headRefName,additions,deletions,files
   ```
2. Extract the linked issue number from `Closes #N` in the PR body.
3. Read the full issue and its comments: `gh issue view <N> --comments` — get the acceptance criterion, references, and any later clarifications.
4. Read PR comments and review comments when present; thread context is part of the review contract, not optional background.
5. Read the epoch roadmap to understand where this fits: `project/epoch_*/roadmap.md`.
6. Read any relevant canonical architecture note in `project/docs/architecture_*.md` when the PR touches public abstractions, ownership boundaries, or recurring patterns.
7. Read `project/components.md` to map the changed files to their expected component scope. Flag if the PR touches files outside that scope without justification.
8. Read only the changed files and the direct dependency files required to review them.
9. Fetch targeted diff context for the touched files. Avoid reading the full PR diff unless the PR is small enough that the broad diff is itself the smallest sufficient context.

## Review lenses

Work through each lens in order. Be specific — name the file and line number when flagging something.

### 1. Mathematical correctness (most important)
- Does the implementation match the mathematical definition cited in the issue?
- Is the right object being computed? (e.g., is the boundary operator actually the transpose of d?)
- Are edge cases handled: empty mesh, boundary cells, degenerate orientations?
- Are the right invariants being tested? Random-input property tests must exist for any new operator.
- Do test names state the invariant explicitly: `test "dd = 0 for ..."`?
- Is floating-point comparison done with explicit epsilon, not exact equality?

### 2. Type safety and comptime
- Is k-form type safety enforced at compile time via `comptime`?
- Can a 2-form be accidentally passed where a 1-form is expected? If so, flag it.
- Are comptime constants used where the value is statically known?
- Are trait/interface requirements checked with `@hasDecl` or `@typeInfo` at comptime?

### 3. Memory and layout
- SoA layout (`std.MultiArrayList`) for any new mesh entity types?
- Explicit allocator passed everywhere — no hidden allocations?
- Are all allocated resources deferred for cleanup (`defer x.deinit(allocator)`)?

### 4. TigerStyle
- `const` by default — any unnecessary `var`?
- All loops bounded — any unbounded iteration?
- Assertions present for preconditions and postconditions (`std.debug.assert`)?
- Comments explain *why*, not what?
- Names are full words, units appended?

### 5. Implementation subtleties
These are the findings the user is least likely to catch:
- Off-by-one errors in index arithmetic
- Sign convention mismatches (orientation of simplices, boundary vs coboundary)
- Numerical stability concerns (division by small quantities, catastrophic cancellation)
- Redundant computation that could be hoisted or cached
- Unnecessary allocations in hot paths
- Code that duplicates logic existing elsewhere in the component
- Weak abstractions introduced under implementation pressure: wrapper nouns,
  builder layers, policy-shaped types, or convenience paths that should have
  been a stronger verb or a direct concept boundary instead
- Cross-component leaks where example or benchmark code is compensating for a
  library seam instead of exercising a coherent public language

### 6. Process
- PR body has `Closes #N`?
- Acceptance criterion from the issue is stated and verified?
- If the PR makes a performance claim, are the benchmark comparisons honest?
- If benchmark methodology changed, was the per-benchmark `version` bumped for the affected rows?
- If benchmark methodology did not change, did the PR preserve base-vs-PR comparability instead of hiding behind a version bump?
- If the claimed perf win comes from a new benchmark with no base counterpart, is there an explicit same-run comparison showing the win?
- Were any non-obvious architectural decisions made? If yes, are they in the decision log?
- If the PR changes the current intended architecture, was the relevant
  canonical note in `project/docs/architecture_*.md` updated?
- Does the PR reveal that a canonical note, design note, or skill definition is
  now stale even if the code change itself is acceptable?
- Labels correct: `type/`, `domain/`, `priority/` all set?
- Does this introduce an interface that conflicts with a known horizon in `project/horizons.md`?
- Do the public names form a coherent abstraction language: stable nouns,
  meaningful verbs, and qualifiers/adjectives in the right place?
- Did this PR miss a chance to remove a weak abstraction that became obvious
  while the seam was already open?
- Are there follow-on issues that should be opened before this merges?
- If review exposes **non-blocking but real** follow-on work that should not bloat the current PR:
  - Check whether it is already tracked
  - If not, recommend opening a **non-milestone** issue with concrete scope and labels
  - When operating autonomously on project hygiene, opening that issue directly is preferred to leaving the work implicit

## Output format

Write a structured review with a clear verdict at the top:

```
## Review: PR #<N> — <title>

**Verdict:** APPROVE | REQUEST CHANGES | BLOCK

### Key findings (implementation detail)
<The most important things the user would miss — specific, with file:line references>

### Mathematical Correctness
<findings>

### Type Safety & Comptime
<findings>

### Memory & Layout
<findings>

### TigerStyle
<findings>

### Implementation Subtleties
<findings>

### Process
<findings>

### Summary
<one paragraph — what's good, what must change before merge>
```

After writing the review, ask the user: "Post this as a GitHub review comment? (y/n)"

If yes:
```sh
gh pr review <number> --comment --body "<review text>"
```
