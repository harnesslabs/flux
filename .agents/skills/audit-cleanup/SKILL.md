---
name: audit-cleanup
description: "Audit the codebase for simplification opportunities: deduplication, unnecessary indirection, API sprawl, shims, over-specialized code, and file-structure cleanup. Use when the repo feels harder to use or change than it should, especially before filing or starting refactor issues."
---

Audit the codebase for simplification pressure. This lens is stricter than style: it asks whether the code has too many concepts, too many public entry points, too many wrappers, or too many files for the capability it delivers.

Usage: /audit-cleanup [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the current library and example surface component-by-component instead of blindly scanning the whole tree.

Before proposing new issues, check the current open GitHub issues for overlapping cleanup/refactor work. Prefer extending or linking existing issues over filing duplicates.

## What to look for

### Public API compression
- Parallel APIs for the same job: a "real" API plus convenience wrappers, aliases, namespace re-exports, or helper overloads
- Public names that leak internal structure instead of presenting one obvious way to do the task
- APIs that ask the caller for information already encoded in types or available from existing values
- Transitional interfaces kept only for backward compatibility; pre-release code should delete them

### Unnecessary indirection
- Wrappers that only rename or forward without enforcing a stronger invariant
- Handle layers, builder layers, or context layers that mostly mirror the underlying operators
- "Generic" adapters whose only current purpose is to hide a direct call path
- Test scaffolding or helper types that exist only because the production API is awkward

### Deduplication and real generic reuse
- Repeated validation logic, slot-management logic, orchestration logic, or serialization logic with only type/degree differences
- Similar examples or modules that should be one parameterized implementation
- Copy-paste APIs that differ only by naming instead of actual semantics
- Genericization opportunities justified by at least two real call sites; do not reward speculative abstraction
- Public APIs that claim to be generic but are implemented as a hardcoded case table over today's supported dimensions, degrees, or variants
- Mathematical or structural definitions that can be expressed once from existing primitives (incidence, recursion, comptime relations) but are instead re-encoded by manual enumeration

### Shim and compatibility removal
- Deprecated pathways, transitional adapters, or "for now" wrappers
- Compatibility code preserving an interface the project no longer wants
- Example-local workarounds that should be either promoted to a real library capability or deleted

### File and module structure
- Files that mix orthogonal responsibilities: core implementation + parsing + example glue + tests for unrelated concerns
- Tiny files whose only job is to bounce into another file
- Modules organized around implementation accidents instead of the user-facing mental model
- Public re-export trees that force users to learn internal folder layout

### Simplicity of use
- Library consumers needing too many steps to do the common thing
- APIs with redundant `withX()` then `x()` patterns when a single demand-driven path would suffice
- Parameter lists carrying compile-time or runtime information that can be derived
- Names that expose mechanics instead of intent
- Generic-looking entry points whose implementations still force maintainers to touch multiple branches for every new degree/dimension case

## Judgment rules

- Favor deletion over abstraction. The best cleanup is often to remove a layer, not generalize it.
- Favor one obvious API. Do not preserve parallel "convenience" and "expert" paths without a current need.
- Favor directness over purity theater. A wrapper that adds no invariant is dead weight.
- Favor generic code only when multiple concrete sites already want the same shape.
- Favor moving capability to the layer where it belongs. Example glue in the library is a smell; library workarounds in examples are also a smell.
- Treat LOC reduction as a proxy, not a goal. Shorter code that weakens invariants is not a win.
- If an API is sold as generic, ask whether the implementation is generic in the same sense. A disguised lookup table is usually not good enough.
- In FEEC/DEC or topology-heavy code, prefer implementations derived from the mathematical object itself (`∂`, incidence, degree recursion, duality relation) over nested switches on `(dimension, degree)`.
- Before accepting a refactor, ask: if we added one more supported degree/dimension, would this code naturally extend, or would we add another case? If another case is needed, call that out.

## Output format

```
## Cleanup Audit: <scope>

Date: YYYY-MM-DD

### High-value simplifications
- <file:line> — <problem>. Why it is unnecessary. Recommended simplification.

### API cleanup
- <file:line> — <extra public surface>. Proposed single obvious API.

### Deduplication / genericization
- <file1:line> and <file2:line> — <repeated pattern>. Shared abstraction or deletion candidate.

### Shims / transitional code
- <file:line> — <compatibility or workaround>. Delete, absorb, or promote.

### Structural cleanup
- <file:line> — <module/file organization problem>. Proposed boundary.

### Recommended issues
| # | Title | Category | Component | Priority |
|---|-------|----------|-----------|----------|
| 1 | ... | api-cleanup | library | high |
```

For each finding, classify it as one of:
- **Fix now** — local deletion, rename, move, or signature simplification with low design risk
- **File/refine issue** — multi-file cleanup, public API change, or a refactor that needs sequencing
- **Keep intentionally** — a layer that looks redundant but currently enforces a real invariant; explain why

Prefer findings that materially simplify the consumer experience or delete maintenance burden. Avoid low-signal nits that belong in `/audit-style`.
