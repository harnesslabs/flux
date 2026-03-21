Check the current epoch and milestone status for the flux project.

Steps:
1. Read `project/epoch_*/roadmap.md` for the most recent epoch — get the current milestone name and its acceptance criterion.
2. Run the following gh commands to get live state:
   ```sh
   gh milestone list --state open
   gh issue list --state open --label "priority/high"
   gh issue list --state open --milestone "<current milestone name>"
   gh issue list --state closed --milestone "<current milestone name>"
   ```
3. Check for tests covering the acceptance criterion: `zig build test`

Report, in order:
1. **Current epoch** — goal and epoch document path
2. **Current milestone** — name, acceptance criterion, issue counts (open / closed)
3. **Blocked issues** — any with `status/blocked` or `status/needs-decision`
4. **High-priority open issues** — list them with their `domain/` and `type/` labels
5. **Acceptance criterion status** — passing or not, based on test output
6. **Recommended next action** — name a specific issue number or implementation step

Label reference for filtering:
- `type/invariant` — math property work (usually on the critical path)
- `status/blocked` — stuck issues that may need a decision
- `priority/high` — must close before milestone is done
