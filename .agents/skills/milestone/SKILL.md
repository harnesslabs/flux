---
name: milestone
description: Plan and create a project milestone with a concrete acceptance criterion and a set of well-scoped capability-level issues.
---

Plan and create a new GitHub milestone for the current epoch.

## Steps

1. Read `project/epoch_N/roadmap.md` for the current epoch to understand existing milestones.
2. Read `project/components.md` to understand scope boundaries for issues.
3. Run `gh milestone list` to see what's already on GitHub.
4. If the user has not stated a milestone goal, ask for it.
5. Propose a single-sentence acceptance criterion — a concrete mathematical or behavioral invariant
   that, when passing in CI, means the milestone is done. Examples:
   - "dd = 0 holds exactly for all k on a uniform 2D triangular mesh with random cochain inputs"
   - "A radiating dipole simulation produces ∇·B = 0 to machine precision at every time step"
   Get explicit confirmation before proceeding.

6. Draft 5–10 GitHub issues. Each issue must:
   - Deliver a **complete capability** (not a single task)
   - Own 3–5 tasks as checkboxes in the issue body
   - At minimum: one test task, one implementation task
   - Require at least one non-obvious design choice
   - Use an imperative verb title stating what capability is delivered:
     "Implement Hodge star ★ for all k-forms with inverse and property tests"
   - Stay within 1–2 component scopes (see `project/components.md`)
   - Have labels assigned:
     - `type/` — what kind of work (`feature`, `invariant`, `test`, etc.)
     - `phase/` — which roadmap phase
     - `domain/` — which codebase area (`topology`, `forms`, `operators`, etc.)
     - `priority/` — `high` if blocking acceptance criterion, else `medium`
   - Include a body with:
     - Description of the capability
     - Acceptance criterion
     - Task checklist (3–5 items)
     - Component scope (from `project/components.md`)
     - References to relevant files, issues, or horizons

   **Sizing check:** If an issue has only 1 task, it's too small — bundle it.
   If it has more than 7 tasks, it's too big — split it.

7. Show the full issue list to the user and get confirmation before creating anything on GitHub.

8. On confirmation:
   a. Create the GitHub milestone:
      ```sh
      gh milestone create --title "<name>" --description "<acceptance criterion>"
      ```
   b. Create each issue with labels:
      ```sh
      gh issue create \
        --title "<title>" \
        --body "<body>" \
        --milestone "<milestone name>" \
        --label "type/feature,phase/N,domain/topology,priority/high"
      ```
   c. Add each issue to the epoch's GitHub Project board:
      ```sh
      gh project item-add <project-number> --owner harnesslabs --url <issue-url>
      ```
      Find the project number with: `gh project list --owner harnesslabs`
   d. Append the milestone to `project/epoch_N/roadmap.md` with its acceptance criterion and milestone URL.

9. Report: milestone URL, issue count, acceptance criterion, and a link to the GitHub Project board.

## Label taxonomy (for reference when assigning)

- type: bug | feature | invariant | perf | docs | refactor | test | ci
- phase: 1 | 2 | 3 | 4 | 5
- domain: topology | forms | operators | io | em | fluid | build
- priority: high | medium | low
