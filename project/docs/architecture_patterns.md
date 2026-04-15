# Architecture Patterns

Date: 2026-04-15

## Purpose

This note records preferred recurring design patterns for flux. The goal is not
to catalog every possible Zig idiom, but to keep the project consistent around a
small number of structural choices that fit the library's design philosophy.

These patterns should be read as defaults. A deviation is allowed, but it
should be justified explicitly.

## Pattern: concept as boundary, not wrapper bait

When a comptime interface check exists, the default design is:

- define the concept
- validate it at the consumer boundary
- pass the conforming type directly

Do **not** automatically introduce a second wrapper type whose only job is to
re-check the concept and forward calls.

In Rust terms, the concept is trait-like. The usual shape should be:

- consumer requires `StepperConcept(Stepper)`
- consumer calls `Stepper.step(...)` directly

not:

- strategy satisfies `StepperConcept`
- `Wrapper(Strategy)` re-validates the same concept
- callers use the wrapper instead of the original type

### Why

A pure forwarding wrapper usually adds no real value:

- no ownership boundary
- no cached state
- no erased representation
- no stronger invariant
- no simplification of the public language

Instead, it creates:

- another noun users must learn
- another layer for tests and examples to route through
- another place for documentation and examples to drift
- more indirection in places that should read directly

### Wrapper admission rule

A wrapper type should exist only if it adds at least one of:

- owned state or lifetime management
- caching or assembled data
- a materially stronger invariant than the underlying type
- a representation change the caller should not see
- a meaningful reduction in user complexity that cannot be achieved by the
  concept alone

If none of those are present, the wrapper is probably dead weight.

## Case study: `TimeStepStrategy` and `TimeStepper`

The current `src/time_stepper.zig` shape likely violates the rule above.

Today:

- `TimeStepStrategy(S)` validates the required `State` and `step` declarations
- `TimeStepper(S)` re-validates the same concept and forwards `State` and `step`

That means the wrapper adds no owned state, no caching, and no stronger
interface. It is just another noun around the same contract.

This then leaks outward:

- tests end up proving that the forwarding wrapper forwards
- examples and helpers start growing indirect “builder to stepper to step”
  structures
- the user-facing language becomes more ceremonial than structural

The stronger default pattern would be:

- keep the concept
- require that concept where a system, runner, or evolution type accepts a
  stepper
- pass the conforming stepper directly

In other words, the interface boundary should be the concept check, not the
wrapper.

## Pattern: test the behavior or the invariant, not the forwarding layer

When a design introduces a forwarding wrapper, it tends to attract tests that
prove plumbing rather than the actual behavior.

Weak:

- “wrapper forwards to inner type”
- “builder returns wrapper”

Stronger:

- the integrator advances state correctly
- the consumer rejects a non-conforming type at comptime
- the PDE invariant holds over repeated steps

This is not an argument against small tests. It is an argument against tests
whose only purpose is to justify an unnecessary layer.

## Pattern: builder types must own real staging complexity

Builder types are allowed, but only when they package real staged setup:

- partially specified configuration
- multi-step assembly with lifetime concerns
- resource ownership that cannot be expressed cleanly with a direct `init`

They should not exist merely to manufacture another thin object whose methods
could have been exposed directly.

If a builder only:

- stores references
- forwards one `initX` call
- returns a thin runtime object

then it is probably an API smell rather than a useful pattern.

## Pattern: examples should reveal the pattern, not compensate for it

Examples are one of the fastest ways to detect a weak standard pattern.

If examples start reading like:

- builder boilerplate
- wrapper boilerplate
- indirection whose only purpose is to satisfy a framework seam

then the framework seam is probably wrong.

The right standard pattern should make examples *shrink*.

## Working rule for future API reviews

Before introducing a new public or semi-public type, ask:

1. Is this a real structural noun?
2. If not, could this be a concept check on an existing noun instead?
3. If it is still a wrapper, what concrete invariant or ownership boundary does
   it add?
4. Would deleting this type make the examples cleaner?

If the answers are weak, do not add the type.

## Current likely follow-up

The `TimeStepStrategy` / `TimeStepper` split should be revisited with this note
in mind. The likely end state is:

- keep the concept
- remove the pure forwarding wrapper
- validate the concept at consumer boundaries such as evolution or system setup

This note does not itself change the implementation. It defines the preferred
pattern so future cleanup and refactors converge toward one consistent shape.

## Pattern: separate model nouns from execution nouns

`System` and `Evolution` are adjacent, but they should not be collapsed by
default.

- `System` is the model/PDE noun
- `Evolution` is the repeated-step execution noun

The right way to simplify this seam is usually not to merge them, but to make
their contract more direct.

That means:

- avoid extra wrapper nouns between the two
- let `Evolution` validate the direct contract it needs from an evolvable system
- keep room for systems that do not evolve and systems that admit multiple
  execution strategies

If a proposed simplification merges `System` and `Evolution`, the burden of
proof is high. It must show that the merge removes real complexity without
hardwiring one execution story into every system.
