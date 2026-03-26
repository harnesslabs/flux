Plan a new epoch for the flux project. This is a back-and-forth planning conversation — not a one-shot command.

## Setup

Before saying anything else:
1. Read all existing `project/epoch_*/roadmap.md` files to understand where we are and what has been done.
2. Read the most recent `project/epoch_*/retrospective.md` — its recommendations feed directly into this planning session.
3. Read `project/initial.md` as the reference architecture.
4. Read `project/horizons.md` to know which future directions must not be precluded by this epoch's work.
5. Read `project/components.md` to understand the current codebase structure.
6. Run `gh api "repos/harnesslabs/flux/milestones?state=all"` and `gh issue list --state open --repo harnesslabs/flux` to see live GitHub state.
7. Present a brief, direct summary of current state (what epoch we're on, what's done, what's open).
8. Ask the user: what should this epoch accomplish? What is the north star?

## During the conversation

Help the user define **2–4 milestones** for the epoch. Each milestone should have enough scope for 5–10 issues (~1–2 weeks of focused AI-assisted work). For each milestone, push until you have:
- A one-sentence goal
- A concrete, testable acceptance criterion (a mathematical invariant or behavioral property)
- A rough list of 5–10 constituent issues to validate scope
- A rough ordering relative to other milestones (dependencies)

### Issue sizing guidelines

Each issue in the rough list should be a **complete capability**, not an atomic task. An issue should:
- Own 3–5 tasks (tests, API design, implementation, edge cases, documentation)
- Require at least one non-obvious design choice
- Produce a PR with ~5+ commits
- Be statable in one sentence: "After this issue, the framework can ___"

**Too small:** "Add boundary check to mesh constructor" — that's a task inside an issue.
**Too big:** "Implement all discrete operators (d, ★, Δ) with composition API" — that's a milestone, not an issue.
**Right-sized:** "Implement Hodge star ★ for all k-forms with inverse and property tests" — clear capability, generic implementation tested at k=0,1,2, multiple tasks, design choices involved.

Push back when:
- A milestone has fewer than 5 natural issues — it may be too small (combine with another)
- A milestone has more than 10 issues — it may be too large (split it)
- An acceptance criterion is vague or not testable
- An issue is a single-task item — bundle it into a larger issue
- Two milestones have an unclear dependency that will cause ordering problems

Ask about the direction beyond the current epoch — architecture decisions made now constrain future work.

## When the plan is agreed

1. Determine the epoch number N from existing `project/` contents.
2. Create the epoch documents:
   - `project/epoch_N/roadmap.md` — structured list of milestones with goals, acceptance criteria, issue lists, and ordering. Include a one-paragraph epoch goal at the top.
   - `project/epoch_N/decision_log.md` — empty template with header only.
   - `project/epoch_N/retrospective.md` — empty template (filled at epoch end via `/retro`).
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

Issues (target 5–10):
- <issue title — one sentence describing the capability>
- <issue title>
...

### M2: <name>
**Goal:** <one sentence>
**Acceptance criterion:** <testable invariant or behavior>
**Depends on:** M1

Issues (target 5–10):
- <issue title>
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

<!-- Written at epoch end via /retro. -->
```
