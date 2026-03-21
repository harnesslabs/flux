Find and work the highest-priority open issue to completion — one issue, one branch, one PR, CI passes, merge.

## Find the issue

1. Read `project/epoch_*/roadmap.md` for the current epoch to understand milestone priorities and acceptance criteria.
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

4. Implement — test first:
   - Write the property test that verifies the acceptance criterion before the implementation
   - Test name must state the invariant: `test "dd = 0 for random 1-forms on triangular mesh"`
   - Fill in the implementation until the test passes

5. Validate locally:
   ```sh
   zig build test   # all tests must pass
   zig fmt src/ build.zig
   zig build        # clean compile
   ```

6. Commit (conventional message, imperative mood):
   ```sh
   git commit -m "feat(topology): implement CSR incidence matrix for 2-simplices"
   ```

## Open the PR

7. Push and open a PR using the project template. Fill every section:
   ```sh
   git push -u origin <branch>
   gh pr create \
     --title "<imperative title>" \
     --body "<filled PR template>" \
     --label "<type/X,phase/N,domain/X,priority/X>"
   ```
   The PR body **must** contain `Closes #<number>` so the issue auto-closes on merge.

8. Report the PR URL. Wait — do not merge until CI passes.

## Merge

9. Check CI: `gh pr checks <number>`
10. When all checks pass:
    ```sh
    gh pr merge <number> --squash --delete-branch
    ```
11. Confirm the issue was auto-closed: `gh issue view <number>`
12. Report: what was merged, what the next recommended issue is.

## Constraints

- Never skip `zig build test` before pushing.
- Never merge with failing CI.
- Never put more than one issue's work in a single PR.
- If during implementation a non-obvious architectural choice is made, run `/decide` before committing.
