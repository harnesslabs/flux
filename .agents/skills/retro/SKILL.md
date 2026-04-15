---
name: retro
description: Review a completed epoch, reconstruct missing architectural decisions, and write the retrospective.
---

Review a completed epoch: reconstruct missing decisions, write the retrospective, and identify process improvements for the next epoch.

Usage: /retro [epoch number]

If no epoch number is given, use the most recent epoch (highest N in `project/epoch_N/`).

## Phase 1: Gather evidence

1. Read the epoch's `roadmap.md` to understand what was planned.
2. Read the epoch's `decision_log.md` to see what was explicitly logged.
3. Read the relevant canonical architecture notes in `project/docs/architecture_*.md` and any targeted design notes in `project/docs/` that the epoch touched.
4. Get the full history of work done in the epoch:
   ```sh
   # All closed issues for each milestone in the epoch
   gh milestone list --state all --json title,number
   # For each milestone:
   gh issue list --state closed --milestone "<milestone>" --json number,title,labels,closedAt
   gh issue list --state open --milestone "<milestone>" --json number,title,labels
   ```
5. Get the PR history:
   ```sh
   gh pr list --state merged --json number,title,mergedAt,additions,deletions,commits
   ```
6. Read the git log for the epoch's timeframe to understand the actual sequence of work:
   ```sh
   git log --oneline --since="<epoch start date>" --until="<epoch end date or now>"
   ```

## Phase 2: Reconstruct decisions

For each merged PR, ask: **was a non-obvious architectural choice made here?**

Signs of an unlogged decision:
- A new type or struct was introduced (implies a modeling choice)
- An existing interface was changed (implies a design tradeoff)
- A dependency between modules was created or removed
- A comment says "we chose X because Y" but the decision log has no entry
- An approach was tried and reverted (visible in git history)

For each unlogged decision found, append it to `project/epoch_N/decision_log.md` with the standard format, tagged `[retroactive]`:

```
## YYYY-MM-DD: <title> [retroactive]

**Decision:** <what was decided>

**Alternatives considered:**
- <alternative and why rejected — reconstruct from PR discussion or code context>

**Rationale:** <best reconstruction of why this choice was made>

**Source:** PR #<N>, commit <sha>
```

Present the reconstructed decisions to the user for review before writing them. The user may correct, add context, or reject entries.

## Phase 3: Write the retrospective

Write `project/epoch_N/retrospective.md` with this structure:

```markdown
# Epoch N Retrospective

Date: YYYY-MM-DD
Duration: <actual time from first commit to last merge>

## What was planned

<Brief summary of the roadmap: milestones, issue counts, acceptance criteria>

## What was delivered

<What actually shipped. Reference PRs and issues by number.>

## Acceptance criteria status

| Milestone | Criterion | Status |
|-----------|-----------|--------|
| M1: ... | ... | ✅ / ❌ / partial |
| M2: ... | ... | ✅ / ❌ / partial |

## What held

<Practices, abstractions, or decisions that worked well and should be carried forward.>

## What didn't hold

<Practices, tools, or assumptions that failed or were never used. Be specific.>

## Process observations

<What was learned about how to work — not what was built, but how it was built. Skill usage, issue sizing, review quality, decision logging, etc.>

## Recommendations for next epoch

<Concrete, actionable changes to skills, CLAUDE.md, issue structure, or workflow. These feed directly into `/epoch` planning.>
```

Present the draft to the user for review. The retrospective is a shared document — the agent writes it, the user edits and approves it.

Alongside the retrospective, propose the canonical doc updates that should be
made to `project/docs/` so the enduring architecture is not trapped only in the
epoch journal.

## Phase 3.5: Close milestones

After the retrospective is written, close all milestones belonging to the epoch on GitHub:

```sh
# For each milestone in the epoch:
gh api repos/:owner/:repo/milestones/<number> -X PATCH -f state=closed
```

Milestones stay open during the epoch so new issues can still be threaded in. Closing them is a retro-time action — it marks the epoch as formally complete.

## Phase 4: Identify improvements

Based on the retrospective, propose specific changes:
- Skill definitions that should be updated (with concrete suggestions)
- CLAUDE.md sections that need revision
- New skills that should be created
- Skills that should be removed or merged
- Changes to issue templates, labels, or GitHub configuration
- Horizon entries that should be added or updated based on what was learned
- Canonical architecture notes that should be created, updated, merged, or retired

These proposals are input to the next `/epoch` session. Do not make the changes directly — present them for discussion.

## Constraints

- Do not fabricate decisions. If you can't determine why a choice was made from the evidence, say so and ask the user.
- The retrospective must be honest. "Everything went great" is not a useful retrospective. Name what didn't work.
- Keep the retrospective concise. It's a working document, not a narrative essay.
- The user has final say on all retrospective content. Present drafts, don't commit final versions without approval.
