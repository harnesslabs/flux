---
name: ideate
description: Explore raw project ideas, classify them, and pressure-test them against the vision, roadmap, decision log, and current architecture.
---

Explore raw ideas with the user and pressure-test them against the project's vision and architecture. This is a back-and-forth conversation — not a one-shot command. The goal is divergent thinking followed by honest filtering, not a plan or commitment.

The argument `$ARGUMENTS` contains the user's raw thoughts, which may be a loose list, a single idea, or a vague direction.

## Setup

Before responding:
1. Read `project/vision.md` to ground yourself in the project's north star.
2. Read `project/horizons.md` to know what future directions are already validated.
3. Read the current epoch's `roadmap.md` (find the latest `project/epoch_N/`) to understand what's already planned.
4. Read the current epoch's `decision_log.md` to know what's already been decided.
5. Read any relevant canonical architecture notes in `project/docs/architecture_*.md` and targeted design notes in `project/docs/` that touch the ideas under discussion.
6. Skim the codebase structure (`src/`) to understand what exists today.

Do NOT dump all of this back at the user. Internalize it silently so your responses are informed.

## Phase 1: Classify

For each idea the user brought, assign it a layer:

| Layer | Meaning | Example |
|-------|---------|---------|
| **Vision** | Changes what flux *is* or *aspires to be* | "support semi-Riemannian metrics" |
| **Architecture** | Changes how the system is structured | "pluggable timesteppers" |
| **Feature** | A concrete capability within the current architecture | "Strang splitting integrator" |
| **Constraint** | A quality or property that cuts across features | "adaptive methods should be trivial" |
| **UX** | How users interact with the system | "structured inputs, no PRM files" |

Present the classification as a concise table. Ask the user if the classification feels right before proceeding.

## Phase 2: Tension-check

For each idea, check it against `project/vision.md`:

- **Aligned**: fits naturally within the stated vision
- **Extends**: goes beyond what the vision says but doesn't contradict it — the vision would need to grow
- **Conflicts**: contradicts something in the vision — one of them would need to change

Be specific about *which* part of the vision it touches. Quote the relevant section.

Flag any ideas that conflict with each other — the user may not have noticed.

## Phase 3: Probe

This is the core of the conversation. For each idea (or cluster of related ideas), push back constructively:

- **Ask for the invariant.** What property should hold if this idea is implemented correctly? If the user can't state one, the idea may not be concrete enough yet.
- **Ask for the cost.** What does this make harder? What complexity does it introduce? What existing simplicity does it sacrifice?
- **Ask for the minimum viable version.** What's the smallest thing that would let you test whether this idea has legs?
- **Connect to existing work.** Does this relate to something already built, planned, or decided? Name the specific file, issue, or decision.
- **Challenge scope.** Is this flux's job, or is it a library/tool that flux should integrate with rather than build?
- **Ask for the abstraction language.** If this idea becomes code, what are the
  nouns, verbs, and adjectives? If the answer is a pile of near-synonyms or
  policy-shaped wrapper types, the idea probably is not cooked enough yet.

Do NOT ask all of these for every idea — use judgment. Some ideas need one question, some need several. Keep the conversation moving.

## Phase 4: Distill

After the back-and-forth has converged (the user will signal this, or you can ask "ready to distill?"), produce a document at `project/ideation/YYYY-MM-DD-<slug>.md` with this structure:

```markdown
# Ideation: <descriptive title>

Date: YYYY-MM-DD

## Context

<1-2 sentences: what prompted this ideation session>

## Ideas explored

### <Idea name>

**Layer:** vision | architecture | feature | constraint | ux
**Vision alignment:** aligned | extends | conflicts
**Summary:** <2-3 sentences capturing the refined version of the idea after discussion>
**Key tension:** <the most important trade-off or open question>
**Verdict:** pursue | park | reject
**Next action:** <if pursue: what concretely happens next — vision amendment, horizons entry, epoch candidate, issue, spike>

### <Idea name>
...

## Connections

<Relationships between ideas that emerged during discussion. Which ones reinforce each other? Which ones are in tension?>

## Open questions

<Questions that weren't resolved and should be revisited>
```

Create the `project/ideation/` directory if it doesn't exist.

If the discussion produces a durable change to the project's architecture
language, ownership model, or recurring design patterns, say so explicitly and
recommend the relevant next write target:
- `project/docs/architecture_*.md` for canonical current design
- `project/horizons.md` for validated but unscheduled future direction
- `project/epoch_N/decision_log.md` for a concrete active decision

## Guidelines

- **Be honest, not encouraging.** If an idea doesn't fit, say so directly. "This is interesting but it's a different project" is a valid and helpful response.
- **Protect simplicity.** Flux's power comes from FEEC/DEC done right, not from being a kitchen-sink framework. Push back on scope creep.
- **Think in types.** Many of these ideas (dimensionful equations, k-form safety) have natural `comptime` expressions in Zig. When an idea maps to a type-level guarantee, say so — that's a strong signal it belongs.
- **No commitments.** This skill produces analysis, not plans. The output is input to `/epoch`, not a substitute for it.
- **Follow the user's energy.** If they're excited about one idea, go deep on it. Don't mechanically march through the list.
- **Keep the canon alive.** If an architecture note or pattern doc is clearly stale relative to the discussion, name that explicitly rather than leaving the drift implicit.
- **Close the loop when obvious.** The skill is still analysis-first, but if the
  user clearly endorses a low-risk documentation or issue-scope update that
  falls directly out of the discussion, make that update instead of stopping one
  step short.
