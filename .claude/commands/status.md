Check the current epoch and milestone status for the flux project.

Steps:
1. Read `project/epoch_*/roadmap.md` for the most recent epoch to understand the current milestone goals and acceptance criteria.
2. Run `gh milestone list --state open` to see active GitHub milestones.
3. Run `gh issue list --state open --milestone <current milestone>` to see open issues.
4. Run `gh issue list --state closed --milestone <current milestone>` to see what's been completed.
5. Check whether any tests exist for the current milestone's acceptance criterion (`zig build test`).

Report:
- Current epoch and its goal
- Current milestone: name, acceptance criterion, deadline if set
- Issues: open vs. closed counts, and which open issues are blocking the acceptance criterion
- Whether the acceptance criterion is currently met (based on test output)
- Recommended next action — be specific: name the next issue or the next implementation step
