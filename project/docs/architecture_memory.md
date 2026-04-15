# Architecture Memory Proposal

Date: 2026-04-15

## Problem

The current project memory model is write-optimized but not read-optimized.

Today, non-obvious decisions are recorded in the current epoch's
`decision_log.md`. That is useful while the work is active, but weak as a
durable reference because it mixes several different kinds of information:

- enduring architectural commitments
- temporary sequencing choices
- issue-local rationale
- experiments that were later superseded

As the project grows, this produces the wrong reading experience:

- too much low-level history to find the durable point
- sparse entries spread across epochs
- repeated context reconstruction when revisiting an area months later

The result is that the decision log is still valuable as a journal, but is not
serving well as the main architecture reference.

## Goal

Keep the discipline of logging non-obvious decisions during active work, while
creating a smaller set of read-optimized documents that preserve the durable
architecture of the project.

The target properties are:

- easy to read after a long context gap
- organized by concept rather than chronology
- explicit about what is current versus superseded
- small enough that agents and humans will actually consult it

## Proposed model

Flux should maintain three distinct documentation layers.

### 1. Journal: epoch decision logs

Location:
- `project/epoch_N/decision_log.md`

Purpose:
- write-optimized record during active work
- capture non-obvious decisions quickly while context is fresh
- preserve alternatives and immediate rationale

What belongs here:
- issue-level design choices
- temporary tradeoffs against the vision
- sequencing choices that matter during the epoch
- local rationale that may later be compressed away

What does not need to survive here forever:
- every transient implementation fork once the durable outcome is known

### 2. Canon: read-optimized architecture notes

Location:
- `project/docs/architecture_*.md`

Purpose:
- durable reference for the current intended design
- concept-oriented, not chronological
- first place to look when re-entering a subsystem

What belongs here:
- stable public nouns and verbs
- ownership boundaries
- cache and invalidation rules
- family distinctions such as DEC vs FEEC
- geometry/metric/topology splits
- benchmark-method policy when it is architectural rather than merely procedural

What does not belong here:
- full local history
- issue-by-issue narration
- every discarded alternative

### 3. Design notes: targeted forward-looking documents

Location:
- `project/docs/<topic>.md`

Purpose:
- focused notes for one architectural thread that is too important to bury in a
  journal entry but not yet broad enough for the canon

Examples:
- wedge / Whitney nonlinear operator direction
- problem setup language
- future geometry/metric API shape

These documents can later either:
- graduate into the architecture canon, or
- remain as scoped references if the topic stays specialized

## Recommended file structure

```text
project/
  docs/
    architecture_overview.md
    architecture_language.md
    architecture_operators.md
    architecture_geometry.md
    architecture_execution.md
    wedge_product_whitney_plan.md
```

The exact split can evolve, but it should stay intentionally small. The point is
not to create a second sprawling doc tree.

## Recommended canonical notes

### `architecture_overview.md`

Purpose:
- short entry point for the durable mental model

Contents:
- what flux is
- the primary abstraction hierarchy
- how to navigate the other architecture notes
- a short “current shape” summary

### `architecture_language.md`

Purpose:
- define the stable public nouns, verbs, and adjectives

Contents:
- `Mesh`, `Geometry`, `Metric`, `System`, `State`, `Evolution` or `Solver`
- criteria for introducing a new public noun
- when qualifiers should be adjectives or comptime parameters
- relationship between examples and public language

This should absorb the durable part of the problem-setup-language ideation.

### `architecture_operators.md`

Purpose:
- define operator-family boundaries and composition rules

Contents:
- DEC versus FEEC roles
- bridge operators
- what belongs in operator contexts
- what should remain family-specific versus family-agnostic

### `architecture_geometry.md`

Purpose:
- define topology, geometry, metric, and motion boundaries

Contents:
- what is purely topological
- what depends on embedding
- what depends on metric
- invalidation rules for moving geometry
- how flat and Riemannian cases relate

### `architecture_execution.md`

Purpose:
- define setup cost, caching, solve/evolution ownership, and benchmark meaning

Contents:
- setup versus hot-path execution
- cache ownership and invalidation
- benchmark selection principles
- how observers, solvers, and integrators fit the model

## Lifecycle

### During an issue

1. Log non-obvious decisions in the current epoch journal.
2. If a targeted topic needs more room, add or update a scoped note in
   `project/docs/`.
3. If the change alters durable architecture, mark the affected canonical note
   for update before the PR is considered complete.

### During milestone or epoch review

1. Distill the journal into the affected canonical notes.
2. Remove or compress entries whose only remaining value is historical.
3. Record superseded or temporary decisions in the retrospective rather than
   keeping them prominent in the canon.

### When a decision is reversed

Do not silently patch the canon and pretend the older design never existed.
Instead:
- update the canonical note to the new current shape
- leave a short “superseded by” pointer in the old scoped note or journal entry
- keep the detailed chronology in the epoch journal or retrospective

This preserves truth without forcing every future reader through the whole
timeline.

## Retention policy

The journal should remain complete within an active epoch, but older journals do
not need to remain the primary reference surface forever.

Recommended policy:

- keep raw epoch logs in place for archival truth
- treat them as append-only historical records, not as the main onboarding path
- periodically distill durable content into the canonical notes
- use retrospectives to summarize what mattered and what did not

This is intentionally conservative: it avoids deleting history while still
making the live architecture easier to read.

## Rules of thumb for what belongs where

Put something in the journal if:
- it was a real decision made under active local context
- the rationale may matter in the next few weeks
- the right answer is still somewhat provisional

Put something in the canon if:
- future work in that area should start from this as the default assumption
- a new contributor or agent would waste time without this context
- the information is organized better by concept than by date

Put something in a scoped design note if:
- it is too large for a journal entry
- it is too specialized or immature for the canon
- multiple future issues are likely to depend on it

## Implications for agent behavior

Agents should stop treating `decision_log.md` as the only architecture memory.
Recommended lookup order for design-sensitive work:

1. `AGENTS.md`
2. `project/components.md`
3. `project/vision.md`
4. relevant canonical architecture note(s) in `project/docs/`
5. current epoch `decision_log.md`
6. targeted scoped note(s) if the topic requires them

That order preserves the current workflow while adding a durable layer between
vision and raw journal history.

## Near-term migration plan

1. Keep writing to `project/epoch_2/decision_log.md` as usual.
2. Promote the durable part of `2026-04-15-problem-setup-language.md` into
   `project/docs/architecture_language.md`. Completed on 2026-04-15.
3. Keep `wedge_product_whitney_plan.md` as a targeted design note.
4. Add `project/docs/architecture_geometry.md` once the topology/geometry/metric
   split is clearer.
5. Update `AGENTS.md` orientation guidance once the first canonical notes exist,
   so the docs are actually consulted.

## Verdict

The project should keep the epoch decision log as a journal, but stop treating
it as the durable architecture reference. The right model is:

- journal for active truth
- canon for current design
- scoped notes for forward-looking or specialized threads

That gives flux a memory structure that is both honest and usable.
