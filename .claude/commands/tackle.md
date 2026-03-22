Find and work the highest-priority open issue to completion — one issue, one branch, one PR, CI passes, merge.

## Find the issue

1. Read `project/epoch_*/roadmap.md` for the current epoch to understand milestone priorities and acceptance criteria.
2. Read `project/horizons.md` — before writing any interface, check that it does not preclude a known future direction.
2. Get live state:
   ```sh
   gh milestone list --state open
   gh issue list --state open --label "priority/high" --json number,title,labels,milestone
   gh issue list --state open --label "status/blocked" --json number,title,labels,milestone
   ```
3. Select the best issue to work next using this priority order:
   - On the current milestone (check which milestone is earliest in the roadmap)
   - `type/invariant` preferred — these are on the critical path for acceptance criteria
   - `priority/high` over `priority/medium`
   - Issues that unblock other `status/blocked` issues
   - Skip `status/blocked` issues (they can't be worked until unblocked)
   - Non-epoch bugs (`type/bug`) may take priority if they are breaking something

4. Present the chosen issue to the user: number, title, acceptance criterion, and a one-sentence rationale for why this one. Ask for confirmation or redirection.

## Work the issue

On confirmation:

1. Read the full issue body: `gh issue view <number>`
2. Create a branch:
   ```sh
   git checkout main && git pull
   git checkout -b <number>-<short-slug>
   ```
   Slug format: lowercase, hyphens, ≤5 words. Example: `42-csr-incidence-matrix`.

3. Read all relevant source files before writing any code. Understand the existing shape before changing it.

4. Implement incrementally — test first, commit often:
   - Write the property test that verifies the acceptance criterion **before** the implementation. Commit: `test(domain): add failing test for <invariant>`.
   - Fill in the implementation in small steps. Each step must compile and pass all existing tests before moving on.
   - **Commit after every coherent unit of progress** — a new type definition, a working constructor, a helper function with its test. A commit is a save state, not a finished feature. Target roughly 50–150 lines per commit; if you are approaching 200+ lines uncommitted, stop and commit what you have.
   - `zig build ci --summary all` between every commit. No exceptions.
   - When the acceptance criterion test passes, that is the final commit for the implementation.

   The rhythm is: **write → test → commit → repeat**. Not: plan everything → write everything → commit once.

## Open the PR

5. Push and open a PR using the project template. Fill every section:
   ```sh
   git push -u origin <branch>
   gh pr create \
     --title "<imperative title>" \
     --body "<filled PR template>" \
     --label "<type/X,phase/N,domain/X,priority/X>"
   ```
   The PR body **must** contain `Closes #<number>` so the issue auto-closes on merge.

6. Report the PR URL and **STOP**. Do not merge. The user must review the PR first.

## After review

Merging is the user's decision. Only merge when the user explicitly says to (e.g., "merge it", "LGTM", "ship it"). Then:

7. Check CI: `gh pr checks <number>` — all checks must pass.
8. Merge:
    ```sh
    gh pr merge <number> --squash --delete-branch
    ```
9. Confirm the issue was auto-closed: `gh issue view <number>`
10. Report: what was merged, what the next recommended issue is.

## Constraints

- **Never merge without explicit user approval.** CI green is necessary but not sufficient. The user reviews every PR.
- Never skip `zig build test` before pushing.
- Never merge with failing CI.
- Never put more than one issue's work in a single PR.
- Never write more than ~200 lines without compiling, testing, and committing. Commits are cheap; debugging a 500-line uncommitted diff is not.
- If during implementation a non-obvious architectural choice is made, run `/decide` before committing.
