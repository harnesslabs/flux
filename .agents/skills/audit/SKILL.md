---
name: audit
description: Run all six audit lenses in parallel and consolidate the findings for this repository or a selected component.
---

Run all six audit lenses in parallel and consolidate the results.

Usage: /audit [component]

If a component is specified (e.g., `operators`), pass it to each sub-audit. Otherwise, audit the full `src/` tree.

## Execution

Spawn six agents in parallel, one for each audit lens:

1. **Safety** — run `/audit-safety` (assertions, bounds, memory, numerical safety)
2. **Style** — run `/audit-style` (naming, dead code, stale references, duplication)
3. **Cleanup** — run `/audit-cleanup` (API sprawl, indirection, shims, deduplication, structure)
4. **Abstraction** — run `/audit-abstraction` (weak nouns and verbs, wrapper layers, policy-shaped APIs, cross-component seams)
5. **Performance** — run `/audit-perf` (layout, allocations, algorithmic efficiency)
6. **Tests** — run `/audit-tests` (coverage, property test quality, redundancy)

Each agent should:
- Read `project/components.md` to understand scope
- Read any relevant canonical architecture note in `project/docs/architecture_*.md` for the scope being audited
- If a component was specified, scope to that component's files and dependencies
- Produce findings in the format defined by its respective skill
- Return its findings to the parent

Across all six lenses, keep one shared abstraction question in view:
- What are the nouns?
- What are the verbs?
- What are the adjectives/qualifiers?

Use that to distinguish true structural problems from isolated local symptoms.

## Consolidation

After all six agents complete, consolidate into a single report:

```
## Full Audit: <scope>

Date: YYYY-MM-DD

### Safety
<summary of critical/warning/note counts>
<top findings>

### Style
<summary>
<top findings>

### Cleanup
<summary>
<top findings>

### Abstraction
<summary>
<top findings>

### Performance
<summary>
<top findings>

### Tests
<summary>
<top findings>

### Recommended issues

| # | Title | Lens | Severity | Component |
|---|-------|------|----------|-----------|
| 1 | ... | abstraction | high | operators/forms |
| 2 | ... | safety | critical | operators |
| 3 | ... | tests | warning | topology |
...
```

Present the consolidated report to the user. Ask which findings they want filed as GitHub issues. For approved findings, create issues with:
- `type/bug` for safety criticals, `type/refactor` for style or cleanup, `type/perf` for performance, `type/test` for test gaps
- Appropriate `domain/` label from the component
- `priority/high` for criticals, `priority/medium` for warnings

If multiple lenses point to the same stale canonical note, targeted design note,
or skill definition, call that out explicitly. Shared process/doc drift is a
real finding, not incidental noise.

## Constraints

- Do not auto-fix anything. The audit produces findings; the user decides what to act on.
- Keep each agent's scope tight — the whole point is that separate lenses need separate context.
- If running on the full codebase, each agent should still work component-by-component to manage context.
