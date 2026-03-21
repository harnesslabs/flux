Log an architectural decision to the current epoch's decision log.

Steps:
1. Identify the current epoch by reading `project/` directory structure.
2. Read the existing `project/epoch_N/decision_log.md` to understand what's already been logged and to avoid duplication.
3. If the user has not provided the full decision context, ask for:
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
```

Today's date is available in the system context. Use it.

Do not editorialize or second-guess the decision. Record it faithfully. The log is a historical record, not a critique.
