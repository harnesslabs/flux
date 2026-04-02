---
name: audit-safety
description: Audit the codebase for safety issues including assertions, bounds, ownership, cleanup, numerical hazards, and invariant enforcement.
---

Audit the codebase for safety issues. Run on `main` between epochs or mid-epoch.

Usage: /audit-safety [component]

If a component is specified (e.g., `operators`), scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

## What to look for

### Assertion coverage
- Every public function should assert its preconditions (`std.debug.assert`)
- Every function that returns computed values should have postcondition assertions where the expected range or property is known
- Loop invariants should be asserted in non-trivial loops
- `comptime` assertions (`comptime { ... }`) for type-level invariants

### Bounds and limits
- All loops must be bounded — flag any loop without a clear termination guarantee
- All queues, buffers, and collections should have explicit capacity limits
- Array indexing should be bounds-checked or provably in-range

### Memory safety
- Every `init`/`create` must have a corresponding `deinit`/`destroy`
- Every allocation must have a `defer` cleanup on the same scope or a documented ownership transfer
- No use-after-free patterns (returning pointers to stack-local data, storing allocator results without tracking)
- No double-free patterns

### Numerical safety
- Division operations should guard against zero divisors (assert or check)
- Subtraction of nearly-equal floats (catastrophic cancellation) should be flagged
- Overflow potential in index arithmetic (especially `usize` operations)

### Invariant enforcement
- Mathematical invariants (dd=0, ∂∂=0, ★★⁻¹=I) must have property tests
- Structural invariants (mesh consistency, orientation agreement) must be asserted

## Output format

```
## Safety Audit: <scope>

Date: YYYY-MM-DD

### Critical (must fix)
- <file:line> — <description of safety issue>

### Warning (should fix)
- <file:line> — <description>

### Note (consider)
- <file:line> — <description>

### Summary
<N critical, M warnings, K notes>
<Recommended issues to create>
```

For each critical or warning finding, propose whether it should be:
- **Auto-fixed now** (obvious fix, no design choice involved)
- **Filed as an issue** (requires design thought or touches multiple files)

Ask the user which findings to act on before creating any issues.
