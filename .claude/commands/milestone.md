Plan and create a new GitHub milestone for the current epoch.

Steps:
1. Read `project/epoch_N/roadmap.md` for the current epoch to understand what milestones already exist and what comes next.
2. If the user has not stated a milestone goal, ask for it.
3. Propose a single-sentence acceptance criterion — a concrete mathematical or behavioral invariant that, when passing in CI, means the milestone is done. Examples:
   - "dd = 0 holds exactly for all k on a uniform 2D triangular mesh with random cochain inputs"
   - "A radiating dipole simulation produces ∇·B = 0 to machine precision at every time step"
   Get explicit confirmation from the user before proceeding.
4. Draft 10–20 GitHub issues that together constitute the milestone work. Each issue should be:
   - Completable in a focused session (2–6 hours)
   - Titled with an imperative verb: "Implement CSR incidence matrix for 2-simplices"
   - Scoped to a single concern — not "implement the Hodge star and test it" (split those)
5. Show the full issue list to the user and get confirmation before creating anything on GitHub.
6. On confirmation:
   a. Create the GitHub milestone: `gh milestone create --title "<name>" --description "<acceptance criterion>"`
   b. Create each issue linked to the milestone: `gh issue create --title "<title>" --body "<body>" --milestone "<milestone name>"`
   c. Append the milestone to `project/epoch_N/roadmap.md` with its acceptance criterion and GitHub milestone link.
7. Report a summary: milestone URL, issue count, and the acceptance criterion.
