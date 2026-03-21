Plan and create a new GitHub milestone for the current epoch.

Steps:
1. Read `project/epoch_N/roadmap.md` for the current epoch to understand existing milestones.
2. Run `gh milestone list` to see what's already on GitHub.
3. If the user has not stated a milestone goal, ask for it.
4. Propose a single-sentence acceptance criterion — a concrete mathematical or behavioral invariant
   that, when passing in CI, means the milestone is done. Examples:
   - "dd = 0 holds exactly for all k on a uniform 2D triangular mesh with random cochain inputs"
   - "A radiating dipole simulation produces ∇·B = 0 to machine precision at every time step"
   Get explicit confirmation before proceeding.
5. Draft 10–20 GitHub issues. Each issue must:
   - Be completable in a focused session (2–6 hours)
   - Use an imperative verb title: "Implement CSR incidence matrix for 2-simplices"
   - Be scoped to one concern (not "implement X and test it" — those are two issues)
   - Have labels assigned:
     - `type/` — what kind of work (`feature`, `invariant`, `test`, etc.)
     - `phase/` — which roadmap phase
     - `domain/` — which codebase area (`topology`, `forms`, `operators`, etc.)
     - `priority/` — `high` if blocking acceptance criterion, else `medium`
   - Include a body with: description, acceptance criterion, and any references
6. Show the full issue list to the user and get confirmation before creating anything on GitHub.
7. On confirmation:
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
        --label "type/feature,phase/1,domain/topology,priority/high"
      ```
   c. Append the milestone to `project/epoch_N/roadmap.md` with its acceptance criterion and milestone URL.
8. Report: milestone URL, issue count, and acceptance criterion.

Label taxonomy (for reference when assigning):
- type: bug | feature | invariant | perf | docs | refactor | test | ci
- phase: 1 | 2 | 3 | 4 | 5
- domain: topology | forms | operators | io | em | fluid | build
- priority: high | medium | low
