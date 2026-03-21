Review a PR for the flux project with mathematical and code correctness scrutiny.

Usage: /review <PR number>

## Gather context

1. Fetch the PR:
   ```sh
   gh pr view <number> --json number,title,body,labels,headRefName,additions,deletions
   gh pr diff <number>
   ```
2. Extract the linked issue number from `Closes #N` in the PR body.
3. Read the full issue: `gh issue view <N>` — get the acceptance criterion and references.
4. Read the epoch roadmap to understand where this fits: `project/epoch_*/roadmap.md`.

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

### 5. Process
- PR body has `Closes #N`?
- Acceptance criterion from the issue is stated in the PR body and verified?
- Were any non-obvious architectural decisions made? If yes, are they in `project/epoch_N/decision_log.md`?
- Labels correct: `type/`, `phase/`, `domain/`, `priority/` all set?

### 6. Molecule checklist
Before approving, verify the PR author answered all downstream questions:
- Does existing documentation need updating (including `project/vision.md`)?
- Does `README.md` need updating to stay consistent with the vision?
- Does `.github/settings.yml` (repo description, topics) need updating?
- Is there a new capability that warrants an end-to-end example?
- Does this change the public API in a way that affects other modules?
- Are there follow-on issues that should be opened before this merges?

A PR that adds new capability without answering these is incomplete, even if CI is green.

## Output format

Write a structured review with a clear verdict at the top:

```
## Review: PR #<N> — <title>

**Verdict:** APPROVE | REQUEST CHANGES | BLOCK

### Mathematical Correctness
<findings — specific, with file:line references>

### Type Safety & Comptime
<findings>

### Memory & Layout
<findings>

### TigerStyle
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
