Audit the codebase for performance issues. Run on `main` between epochs or mid-epoch.

Usage: /audit-perf [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

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

Ask the user which findings warrant issues.
