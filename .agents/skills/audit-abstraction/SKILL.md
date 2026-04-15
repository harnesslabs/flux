---
name: audit-abstraction
description: Audit the codebase for weak abstractions, policy-shaped APIs, redundant wrappers, and cross-component seams that should be simplified by deletion, consolidation, or stronger structural nouns and verbs.
---

Audit the codebase for abstraction quality. This lens is broader than cleanup:
it looks for places where the public or internal design has drifted away from a
coherent domain language and should be simplified, generalized, or deleted.

Usage: /audit-abstraction [component]

If a component is specified, start there per `project/components.md`, but do
not stop at the boundary if the abstraction problem is created by duplicated or
conflicting seams across components. Cross-component findings are allowed when
they materially simplify the tree.

Before proposing new issues, check the current open GitHub issues for overlap.
Prefer extending or linking existing cleanup or refactor issues over filing
duplicates.

Before diving into code, read the relevant canonical architecture notes in
`project/docs/architecture_*.md` and any targeted design note in `project/docs/`
for the area being audited. Audit against the intended language, not only the
current implementation.

## First principles

- Favor deletion over addition. The default fix for a weak abstraction is to
  remove a layer or collapse duplicate concepts, not to introduce a wrapper.
- Favor one obvious way. Parallel APIs, transitional convenience paths, and
  duplicate public nouns are presumed guilty until justified.
- Favor structural nouns and strong verbs. Policy-shaped wrappers and
  case-specific names are usually evidence that the real abstraction has not
  been named yet.
- Favor generalization that removes code. Do not reward “generic” refactors
  that add indirection, forwarding, or more surface area than they delete.
- Strike while the iron is hot. If you are already in a seam and can see a
  stronger local API that simplifies the area you are touching, surface that
  opportunity immediately instead of paving over it for later.
- Keep shared memory alive. If the audit reveals a durable reusable pattern,
  stale canonical note, or missing architecture note, say so explicitly and
  recommend the write target in `project/docs/`.

## What to look for

### Weak public language
- Public APIs whose nouns, verbs, and adjectives do not form a small coherent
  domain language
- Near-synonym public types that should be one noun with adjectives or one noun
  plus a stronger verb
- Builder or context layers that mirror underlying operations without imposing
  a stronger invariant
- API surfaces that expose implementation accidents rather than mathematical or
  computational structure

### Additive cleanup smells
- A “cleanup” that added a wrapper instead of removing a weak seam
- Compatibility aliases, shims, or helper overloads in pre-release code
- Convenience and expert APIs that both express the same capability
- Public names added only to avoid deleting an older one

### Cross-component abstraction drift
- Duplication created at component boundaries because the boundary is drawn in
  the wrong place
- Similar ownership models expressed differently across neighboring modules
- Examples or benchmarks carrying concepts that should be library-level nouns
- Library abstractions that exist only to compensate for awkward example or
  benchmark entry points

### False genericity
- APIs presented as general but implemented by manual case tables over today's
  dimensions, degrees, families, or variants
- Generic wrappers with only one meaningful consumer
- Abstractions that would require another branch, noun, or file family as soon
  as one more PDE, degree, or mesh variant is added

### Deletion candidates
- Deprecated, fallback, legacy, or transitional code
- Thin forwarding helpers with no new invariant
- Stubs, placeholder surfaces, AI-generated filler, or comments narrating code
  churn instead of helping a new reader
- Tests that only exist because the API is awkward rather than because they
  protect a unique invariant

## Pressure tests

For each suspected abstraction problem, ask:

- What should cease to exist if we did this right?
- What is the stable noun here?
- What is the real verb?
- Which qualifiers are adjectives rather than reasons to mint a new noun?
- If we add one more supported case tomorrow, would the design extend naturally
  or would we add another wrapper, switch branch, or sibling type?
- Is this actually flux's job, or is the code compensating for a missing lower
  layer or an overly opinionated upper layer?

If the likely fix adds code, justify why deletion or consolidation would fail.
If that justification is weak, the finding should recommend removal instead.

When the right pattern depends on broader Zig practice and local repo context is
not enough, do targeted external research before making a strong recommendation.
Prefer:
- Zig language and standard library documentation
- established Zig projects whose style is relevant to the pattern under review
- a small number of primary sources rather than broad blog-search churn

Bring the result back into flux's own language. Do not cargo-cult a community
pattern if it conflicts with the project's architecture notes.

## Output format

```
## Abstraction Audit: <scope>

Date: YYYY-MM-DD

### High-confidence deletions
- <file:line> — <layer or path>. Why it should cease to exist. Proposed simpler shape.

### Consolidations
- <file1:line> and <file2:line> — <duplicated concept>. Proposed single noun/verb boundary.

### Cross-component seams
- <component(s)> / <file:line> — <boundary problem>. Why the current split creates duplication or weak APIs.

### Keep intentionally
- <file:line> — <layer that looks redundant>. What invariant it actually enforces.

### Shared-memory updates
- <doc or skill path> — <what is stale or missing>. Recommended update.

### Recommended issues
| # | Title | Category | Component | Priority |
|---|-------|----------|-----------|----------|
| 1 | ... | abstraction-cleanup | operators/forms | high |
```

For each finding, classify it as one of:
- **Fix now** — local deletion, consolidation, rename, or seam tightening with
  low design risk
- **File/refine issue** — multi-file or public abstraction change that needs
  sequencing
- **Keep intentionally** — the layer enforces a real invariant and should stay

Prefer a few high-signal findings over a long list of nits.

When the audit converges on an obvious low-risk shared-memory or issue update,
perform it by default instead of stopping at recommendation only. Typical
examples:
- extend an existing issue whose scope clearly should absorb the finding
- update a canonical architecture note to capture the conclusion
- patch a stale skill definition that the audit directly exposed

Only stop at recommendation when the write target is genuinely ambiguous or the
change would commit the project to a direction the user has not yet endorsed.

## Relationship to other audit lenses

- Use `/audit-cleanup` for broad simplification, deduplication, and file/module
  cleanup.
- Use this skill when the main question is whether the *abstraction language*
  itself is wrong, too additive, too wrapper-heavy, or split at the wrong seam.
- If a finding is mostly about naming or comments, it belongs in
  `/audit-style`.
- If a finding is mostly about unused code or performance, route it to the
  appropriate lens and keep this audit focused on structure.
