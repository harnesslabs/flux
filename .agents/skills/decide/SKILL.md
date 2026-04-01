---
name: decide
description: Record a non-obvious architectural decision in the current epoch decision log with alternatives and rationale.
---

Log an architectural decision to the current epoch's decision log.

This skill is primarily **agent-invoked during `/tackle`** when a non-obvious design choice is made. It can also be invoked by the user at any time. The agent should call this whenever it encounters a fork in the road during implementation — struct layout, ownership model, comptime parameter choices, interface shape, algorithm selection.

**Trigger signals** (the agent should recognize these during `/tackle`):
- Choosing between two reasonable struct layouts
- Deciding who owns an allocation (caller vs callee)
- Picking a comptime parameter strategy over a runtime one (or vice versa)
- Selecting an algorithm when multiple valid options exist
- Changing an existing interface to accommodate new requirements
- Deviating from a pattern established elsewhere in the codebase

## Steps

1. Identify the current epoch by reading `project/` directory structure.
2. Read the existing `project/epoch_N/decision_log.md` to understand what's already been logged and to avoid duplication.
3. If the context is not already clear (e.g., the user invoked this directly without prior discussion), ask for:
   - The decision made (what was chosen)
   - Alternatives that were considered (at least one)
   - The rationale — why this option over the alternatives
4. Append the decision to `project/epoch_N/decision_log.md` in this format:

```
## YYYY-MM-DD: <short imperative title>

**Decision:** <what was decided, one or two sentences>

**Alternatives considered:**
- <alternative 1 and why it was rejected>
- <alternative 2 and why it was rejected>

**Rationale:** <the reasoning — constraints, trade-offs, future implications>

**Source:** <PR #N or issue #N if applicable>
```

Today's date is available in the system context. Use it.

Do not editorialize or second-guess the decision. Record it faithfully. The log is a historical record, not a critique.

When invoked during `/tackle`, commit the decision log update alongside the code it pertains to — don't make a separate commit just for the log entry.
