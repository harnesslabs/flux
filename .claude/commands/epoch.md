Plan a new epoch for the flux project. This is a back-and-forth planning conversation — not a one-shot command.

## Setup

Before saying anything else:
1. Read all existing `project/epoch_*/roadmap.md` files to understand where we are and what has been done.
2. Read `project/initial.md` as the reference architecture.
3. Read `project/horizons.md` to know which future directions must not be precluded by this epoch's work.
3. Run `gh api "repos/harnesslabs/flux/milestones?state=all"` and `gh issue list --state open --repo harnesslabs/flux` to see live GitHub state.
4. Present a brief, direct summary of current state (what epoch we're on, what's done, what's open).
5. Ask the user: what should this epoch accomplish? What is the north star?

## During the conversation

Help the user define **3–5 milestones** for the epoch. Each milestone should have enough scope for 10–20 issues (~1 week of focused AI-assisted work). For each milestone, push until you have:
- A one-sentence goal
- A concrete, testable acceptance criterion (a mathematical invariant or behavioral property)
- A rough list of 10–15 constituent issues to validate scope
- A rough ordering relative to other milestones (dependencies)

Push back when:
- A milestone is too small (fewer than ~10 natural issues) — consolidate it
- A milestone is too large to close in ~1 week — split it
- An acceptance criterion is vague or not testable
- Two milestones have an unclear dependency that will cause ordering problems

Ask about the direction beyond the current phase — architecture decisions made now constrain future phases.

## When the plan is agreed

1. Determine the epoch number N from existing `project/` contents.
2. Create the epoch documents:
   - `project/epoch_N/roadmap.md` — structured list of milestones with goals, acceptance criteria, issue lists, and ordering. Include a one-paragraph epoch goal at the top.
   - `project/epoch_N/decision_log.md` — empty template with header only.
   - `project/epoch_N/retrospective.md` — empty template (filled at epoch end).
3. Create a GitHub issue to track the planning artifact:
   ```sh
   gh issue create \
     --repo harnesslabs/flux \
     --title "docs: plan epoch N — <title>" \
     --body "..." \
     --label "type/docs,priority/high,domain/build"
   ```
   The body should summarize the epoch goal, list the milestones, and state the epoch acceptance criterion.
4. Create a branch for that issue:
   ```sh
   git checkout -b {issue-number}-plan-epoch-N
   ```
5. Commit the epoch documents on that branch:
   ```sh
   git add project/epoch_N/
   git commit -m "docs: add epoch N roadmap and planning documents"
   ```
6. Open a draft PR:
   ```sh
   gh pr create --draft \
     --repo harnesslabs/flux \
     --title "docs: epoch N — <title>" \
     --body "..."
   ```
7. Report the epoch document paths, issue URL, and draft PR URL.
8. Tell the user to run `/milestone` for each milestone when ready to create GitHub Milestones and issues.

## Templates

### roadmap.md
```markdown
# Epoch N: <title>

<One paragraph: what this epoch accomplishes and why it matters.>

## Milestones

### M1: <name>
**Goal:** <one sentence>
**Acceptance criterion:** <testable invariant or behavior>
**Depends on:** none

Issues (target 12–15):
- <issue 1>
- <issue 2>
...

### M2: <name>
**Goal:** <one sentence>
**Acceptance criterion:** <testable invariant or behavior>
**Depends on:** M1

Issues (target 12–15):
- <issue 1>
...
```

### decision_log.md
```markdown
# Epoch N Decision Log

<!-- Append entries with /decide. Format: ## YYYY-MM-DD: <title> -->
```

### retrospective.md
```markdown
# Epoch N Retrospective

<!-- Written at epoch end. What held? What didn't? What would we change? -->
```
