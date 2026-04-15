---
name: audit-perf
description: Audit the codebase for performance issues such as memory layout, allocation patterns, hot loops, and sparse operator inefficiencies.
---

Audit the codebase for performance issues. Run on `main` between epochs or mid-epoch.

Usage: /audit-perf [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

Read the relevant canonical architecture note in `project/docs/architecture_*.md`
and benchmark-policy guidance in `AGENTS.md` before recommending structural perf
changes. Performance advice should reinforce the intended execution model.

## What to look for

### Memory layout
- SoA layout (`std.MultiArrayList`) for entity data — flag AoS patterns
- Cache-friendly access patterns in hot loops (sequential access, not random)
- Unnecessary copies where slices or references would suffice

### Allocation patterns
- Allocations inside loops that could be hoisted (allocate once, reuse)
- Temporary allocations that could use a fixed-buffer or arena allocator
- Growing arrays that could be pre-sized with `ensureTotalCapacity`

### Algorithmic efficiency
- O(n²) patterns where O(n log n) or O(n) alternatives exist
- Redundant computation (same value computed multiple times in a loop)
- Unnecessary sorting, searching, or iteration

### Sparse matrix operations
- CSR/CSC access patterns that cause cache misses
- Dense operations on sparse data
- Matrix assembly that could be done in-place

### SIMD and vectorization readiness
- Data layouts that prevent auto-vectorization (AoS, scattered access)
- Branch-heavy inner loops that block vectorization
- Operations on small fixed-size arrays that could use `@Vector`

### Hot path identification
- Which functions are likely hot paths (called per-timestep, per-cell, per-edge)?
- Are hot paths free of allocations and assertions that could be `@setRuntimeSafety(false)` in release?

### Abstraction/performance tension
- Generic-looking surfaces whose adjective space (`dimension`, `degree`, scalar type, layout mode) explodes into manual runtime branching
- Wrapper layers that hide data movement, allocation, or copies behind tidy API names
- Policy-specific nouns that prevent reuse of one optimized hot path across several variants

## Output format

```
## Performance Audit: <scope>

Date: YYYY-MM-DD

### Hot paths
<List the functions identified as likely hot paths and their current characteristics>

### Critical (clear performance bug)
- <file:line> — <description, estimated impact>

### Opportunity (measurable improvement possible)
- <file:line> — <description, what to change>

### Future (relevant when scale increases)
- <file:line> — <description, when this matters>

### Summary
<overall assessment, recommended issues to create>
```

All findings should be filed as issues, not auto-fixed. Performance changes require benchmarking to verify improvement, which means they need their own PR with before/after measurements.

When a performance issue is really an abstraction issue, say so explicitly:
- wrong noun: too many concrete wrapper types block one optimized path
- wrong verb: work is split across multiple steps instead of one fused operation
- wrong adjective: a qualifier that should be `comptime` or layout-level is handled as dynamic casework

Benchmark hygiene is part of the audit:
- Flag apples-to-oranges comparisons where a benchmark method changed without a per-benchmark `version` bump.
- Flag unnecessary whole-suite invalidation when only one benchmark row changed methodology.
- Prefer durable evidence: unchanged benchmark rows for base-vs-PR comparisons, plus same-run comparison benchmarks when a new benchmark is introduced.

Ask the user which findings warrant issues.

If the hot paths the audit identifies contradict the project's current setup or
execution language, say so explicitly and name the canonical note that should be
updated.
