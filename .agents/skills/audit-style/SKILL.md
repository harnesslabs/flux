---
name: audit-style
description: Audit the codebase for naming, hygiene, comment quality, duplication, formatting, and structural consistency issues.
---

Audit the codebase for style and hygiene issues. Run on `main` between epochs or mid-epoch.

Usage: /audit-style [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

## What to look for

### Naming
- Full words, no abbreviations: `exterior_derivative` not `ext_deriv`
- Units and qualifiers appended in descending significance
- Consistent naming patterns across similar constructs
- Public names should read like a coherent domain language: stable nouns,
  clear verbs, and qualifiers/adjectives that refine rather than fragment the model

### Dead code and stale references
- Unreachable code paths
- Commented-out code blocks
- Inline comments referencing issue numbers, PR numbers, or other artifacts that may be stale
- TODO/FIXME/HACK comments that should be filed as issues or resolved
- Unused imports, variables, or functions

### Comment quality
- Comments should explain *why*, not *what*
- Mathematical content should use Unicode symbols (∂, ∇, ∈, ★, ∧, ⟨⟩, Ω)
- Citations should include author, title, and equation number
- No fabricated URLs

### Code duplication
- Similar logic repeated in multiple places that could share a common implementation
- Copy-paste patterns with minor variations

### Structural consistency
- Similar constructs should follow similar patterns (e.g., all operator types should have the same method shape)
- Public API should be consistent in style (all snake_case, consistent parameter ordering)
- Nested `switch` / `if` ladders over semantic axes like `(dimension, degree)` that may indicate a missing derived definition or comptime formulation rather than a true need for casework
- Generic-looking APIs whose implementation style contradicts the abstraction they present
- Families of type names that differ only by policy adjectives and may want one stronger noun instead

### Formatting
- `zig fmt` compliance (should be caught by CI, but verify)
- Consistent file organization (imports, types, functions, tests)

## Pressure tests

Before concluding that a style issue is only cosmetic, ask:

- Is this "generic" function actually generic, or is it a manually enumerated case table?
- Could the implementation be written once in terms of the native structure (incidence, recursion, comptime relation) instead of today’s supported cases?
- Does this control flow naturally extend when the next degree/dimension/variant is added, or does it require editing multiple branches?
- Are these words the right nouns, verbs, and adjectives for the abstraction, or are the names exposing implementation accidents?

## Output format

```
## Style Audit: <scope>

Date: YYYY-MM-DD

### Stale references
- <file:line> — <what's stale and why>

### Dead code
- <file:line> — <what's dead>

### Naming issues
- <file:line> — <current name> → <suggested name>

### Duplication
- <file1:line> and <file2:line> — <description of duplicated logic>

### Comment issues
- <file:line> — <what's wrong with the comment>

### Summary
<overall assessment, recommended issues to create>
```

For each finding, propose whether it should be:
- **Auto-fixed now** (rename, remove dead code, fix a comment — no design choice)
- **Filed as an issue** (duplication that needs a shared abstraction, structural inconsistency)

Ask the user which findings to act on before making changes or creating issues.
