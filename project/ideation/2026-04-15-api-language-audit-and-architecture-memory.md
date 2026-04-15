# Ideation: API Language, Abstraction Audit, and Architecture Memory

Date: 2026-04-15

## Context

This ideation session started from a backlog of notes collected during recent
example and API cleanup work. The discussion focused on which notes still
describe real design pressure, which were already resolved by newer
abstractions, and which should drive the next project artifacts.

## Ideas explored

### Mesh-owned operators

**Layer:** architecture
**Vision alignment:** conflicts
**Summary:** The original note suggested that operators should become methods on
the mesh or otherwise be mesh-owned because operators feel intrinsic to the
mesh. After review, this was judged stale. The current direction with explicit
operator contexts is stronger because it keeps topology ownership separate from
assembled, family-specific, and problem-specific operator state.
**Key tension:** A mesh-centric user experience is desirable, but mesh-owned
numerical state would blur topology, geometry, metric dependence, and assembly
policy in exactly the place the architecture needs explicit separation.
**Verdict:** reject
**Next action:** None. Preserve the current operator/context ownership model.

### FEEC and DEC as explicit operator families

**Layer:** architecture
**Vision alignment:** extends
**Summary:** FEEC and DEC should remain explicit public mathematical families,
with room for mixed formulations where strong-form DEC is the right tool in one
part of a problem and weak-form FEEC is the right tool elsewhere. The
distinction should remain semantic rather than historical: both families can be
present, but not as duplicate legacy surfaces for the same operation.
**Key tension:** The project wants both mathematical honesty and one obvious way
to do things. That means mixed FEEC/DEC composition is acceptable only when the
underlying semantics genuinely differ.
**Verdict:** pursue
**Next action:** Keep pressure-testing milestone work against the question: does
mixed user code read like one coherent domain language rather than glue between
two sub-libraries?

### Explicit metric and geometry semantics

**Layer:** constraint
**Vision alignment:** aligned
**Summary:** The discussion converged on explicit structural nouns such as
`Metric(.euclidean)` and `Metric(.riemannian, g)` rather than vague mode flags.
At the same time, this likely points to a broader split between topology,
geometry/embedding, and metric-dependent structure. Euclidean geometry is a
special case, but the architecture should not flatten that distinction too
early.
**Key tension:** A clean `Metric` noun is attractive, but moving meshes and ALE
support may require `Metric` to sit inside a larger `Geometry` or equivalent
object rather than carry all geometric meaning itself.
**Verdict:** pursue
**Next action:** Use the problem-setup design exercise to determine the right
top-level nouns before committing to an API.

### Problem-setup-first API design

**Layer:** ux
**Vision alignment:** aligned
**Summary:** The strongest insight from the conversation is that the project has
sometimes been designing APIs backward from implementation and examples instead
of first writing the ideal simulation-setup surface. The next design artifact
should therefore be a note with 2-3 idealized setup sketches that force the
project to name its structural nouns, verbs, compile-time parameters, cache
boundaries, and intended hot paths.
**Key tension:** This is not merely documentation work. If treated as docs too
early, it will describe the current API instead of revealing the better one.
**Verdict:** pursue
**Next action:** Create a problem-setup design note with idealized API sketches
for fixed-mesh Maxwell, a metric-aware elliptic/parabolic problem, and one
moving-geometry or ALE-flavored toy problem.

### Benchmark design from user hot paths

**Layer:** constraint
**Vision alignment:** aligned
**Summary:** Benchmark priorities should come from intended user workflows and
dominant execution paths, not from whichever functions are easiest to time.
This ties benchmark design to the same setup-language exercise: the code path
that matters most to the user should be visible from the problem setup and
execution model.
**Key tension:** Optimization is not the top priority, but benchmark choices
shape engineering attention. Measuring implementation accidents instead of
product-critical paths would quietly distort the library.
**Verdict:** pursue
**Next action:** When the setup sketches exist, derive benchmark targets from
their hot paths and setup costs.

### Abstraction-audit skill focused on deletion and consolidation

**Layer:** architecture
**Vision alignment:** extends
**Summary:** The existing audit lenses are useful, but the project needs a
dedicated abstraction-focused audit that maps public nouns, verbs, and
ownership boundaries, then defaults toward simplification by deletion rather
than additive wrappers. This should explicitly fight the common LLM failure mode
of cleaning up by introducing another layer instead of removing the weak one.
**Key tension:** Cleanup often crosses component boundaries, while the current
scope discipline intentionally narrows local context. The skill must preserve
focus without becoming blind to duplication or bad seams created at component
boundaries.
**Verdict:** pursue
**Next action:** Design an abstraction-audit skill and update broader agent
guidance so in-flight implementation work also prefers removal, consolidation,
and immediate local API strengthening when the seam is already open.

### Redundant tests

**Layer:** constraint
**Vision alignment:** aligned
**Summary:** The concern is not with the existence of many tests, but with truly
redundant tests that restate weaker guarantees and increase maintenance cost.
The distinction between mathematical proof obligations, API behavior checks, and
optimized-path regressions remains important.
**Key tension:** Aggressive pruning can accidentally delete useful coverage if
the suite is not first classified by what invariant each test uniquely protects.
**Verdict:** park
**Next action:** Revisit after the abstraction-audit and setup-language work,
when the intended surfaces and invariants are clearer.

### Read-optimized architecture memory

**Layer:** architecture
**Vision alignment:** conflicts
**Summary:** The current decision log is useful as a write-optimized journal but
weak as durable project memory. The discussion converged on keeping the raw log
during active work while periodically distilling durable outcomes into a smaller
read-optimized architecture record.
**Key tension:** The current project doctrine says every non-obvious decision is
logged. The improvement is not to stop logging, but to stop treating the raw
journal as the primary long-term reference.
**Verdict:** pursue
**Next action:** After the setup-language artifact and abstraction-audit design,
draft a proposal for a distilled architecture record and the retention/compression
policy for epoch logs.

## Connections

The discussion converged on one primary dependency chain:

1. Write the ideal problem-setup sketches first.
2. Use those sketches to define the abstraction-audit lens and its deletion-first
   philosophy.
3. Use both of those to decide what durable architectural knowledge should be
   preserved in a read-optimized form.

Metric/geometry semantics, FEEC/DEC composition, and benchmark priorities all
feed into the setup-language exercise. The abstraction-audit work then turns the
resulting language into a cleanup and simplification discipline. The architecture
memory proposal comes last because it should preserve the refined language, not
the raw exploratory churn.

## Open questions

- Is `Metric` the right top-level noun, or should it be one field inside a
  broader `Geometry` object?
- For mixed FEEC/DEC problems, what are the stable user-facing nouns and verbs
  that keep the composition coherent instead of feeling like interop glue?
- How should abstraction-focused cleanup intentionally cross component
  boundaries without collapsing back into “scan the whole repo” behavior?
- Which agent instructions should change globally so deletion and consolidation
  are default cleanup moves during ordinary implementation, not only during a
  dedicated audit?
