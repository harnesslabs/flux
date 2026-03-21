Plan a new epoch for the flux project. This is a back-and-forth planning conversation — not a one-shot command.

## Setup

Before saying anything else:
1. Read all existing `project/epoch_*/roadmap.md` files to understand where we are and what has been done.
2. Read `project/initial.md` as the reference architecture.
3. Run `gh milestone list --state all` and `gh issue list --state open` to see live GitHub state.
4. Present a brief, direct summary of current state (what epoch we're on, what's done, what's open).
5. Ask the user: what should this epoch accomplish? What is the north star?

## During the conversation

Help the user define **5–10 milestones** for the epoch. For each milestone, push until you have:
- A one-sentence goal
- A concrete, testable acceptance criterion (a mathematical invariant or behavioral property)
- A rough ordering relative to other milestones (dependencies)

Push back when:
- A milestone is too large to close in ~1 week of focused AI-assisted work
- An acceptance criterion is vague or not testable
- Two milestones have an unclear dependency that will cause ordering problems

Ask about the direction beyond the current phase — architecture decisions made now constrain future phases.

## When the plan is agreed

1. Determine the epoch number N from existing `project/` contents.
2. Create the epoch documents:
   - `project/epoch_N/roadmap.md` — structured list of milestones with goals, acceptance criteria, and ordering. Include a one-paragraph epoch goal at the top.
   - `project/epoch_N/decision_log.md` — empty template with header only.
   - `project/epoch_N/retrospective.md` — empty template (filled at epoch end).
3. Create a GitHub Project for this epoch:
   ```sh
   gh project create --owner harnesslabs --title "Epoch N: <title>"
   ```
4. Report the epoch document paths and GitHub Project URL.
5. Tell the user to run `/milestone` for each milestone when ready to create GitHub Milestones and issues.

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

### M2: <name>
**Goal:** <one sentence>
**Acceptance criterion:** <testable invariant or behavior>
**Depends on:** M1
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
