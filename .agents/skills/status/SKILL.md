---
name: status
description: Summarize the current epoch, milestone, open work, blockers, and acceptance-criterion status from docs, GitHub state, and tests.
---

Check the current epoch and milestone status for the flux project.

## Steps

1. Read `project/epoch_*/roadmap.md` for the most recent epoch — get the current milestone name and its acceptance criterion.
2. Read relevant canonical architecture notes in `project/docs/architecture_*.md` if current work is clearly clustered around one architectural thread.
3. Read `project/components.md` for component context.
4. Run the following gh commands to get live state:
   ```sh
   gh milestone list --state open
   gh issue list --state open --label "priority/high"
   gh issue list --state open --milestone "<current milestone name>"
   gh issue list --state closed --milestone "<current milestone name>"
   ```
5. Check for tests covering the acceptance criterion: `zig build test --summary all`

Report, in order:
1. **Current epoch** — goal and epoch document path
2. **Current milestone** — name, acceptance criterion, issue counts (open / closed)
3. **In-progress work** — any open PRs or draft PRs
4. **Blocked issues** — any with `status/blocked` or `status/needs-decision`
5. **High-priority open issues** — list them with their `domain/` and `type/` labels
6. **Acceptance criterion status** — passing or not, based on test output
7. **Decision log health** — are there entries since last check? If the milestone has had PRs merged but no decisions logged, flag it.
8. **Architecture doc health** — do the canonical notes in `project/docs/` look current for the active work, or is there obvious drift?
9. **Recommended next action** — name a specific issue number or suggest `/tackle`

Label reference for filtering:
- `type/invariant` — math property work (usually on the critical path)
- `status/blocked` — stuck issues that may need a decision
- `priority/high` — must close before milestone is done
