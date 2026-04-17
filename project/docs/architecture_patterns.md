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

## Case study: execution stepper cleanup

Flux previously carried a split between a stepper concept and a pure forwarding
wrapper. That shape was removed because it violated the rule above.

The old shape:

- re-validated the same contract in multiple layers
- taught wrapper-oriented tests
- leaked adapter and builder boilerplate into examples

The current preferred pattern is:

- keep one direct method concept at the execution boundary
- let `Evolution` validate the chosen method directly
- let `Evolution` own `dt`, time, and run bookkeeping
- route method-specific configuration through `Evolution.Config`

In other words, the interface boundary is the concept check plus the execution
owner, not a wrapper around the concept.

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

## Pattern: type factories remain nouns; initialization lives on the returned type

Many flux nouns are legitimately parameterized at comptime. That means a
type-level factory such as:

- `Mesh(2, 2)`
- `Cochain(MeshType, k, Duality)`
- `Evolution(StateType, MethodType, AuxType)`

may still be the right implementation tool and the right public noun.

The preferred construction pattern is:

- keep the noun as the type factory
- put staged configuration and initialization on the returned type
- avoid module-level `init(...)` helpers that make the namespace pretend to be
  the noun

### Why

This matches the mesh pattern cleanly:

- `Mesh(...).plane(...)`
- `Evolution(...).config().dt(...).init(...)`

Both read as:

- name the structural noun
- choose construction/configuration on that noun
- produce the runtime value

The anti-pattern is:

- `module.init(...)`

because it weakens the package language and makes the namespace, rather than
the noun, look like the constructor surface.

Examples should also avoid redundant type aliases such as:

- `const EvolutionType = Evolution(...);`
- `var evolution = try EvolutionType.config()...init(...);`

unless the type name is genuinely needed in type position.

## Pattern: avoid fake runtime nouns when a method is really type-level

If an execution method has no independent runtime identity, keep it as a type.

Weak:

- runtime `Stepper` object stores `state`
- runtime `Stepper` object stores `dt`
- `Evolution` also stores `state` or run counters

Stronger:

- `Method` is a comptime type with `advance(...)`
- `Evolution` owns `dt`, time, and run bookkeeping
- method-specific settings flow through `Evolution.Config`

This keeps the execution language flat and avoids double ownership of the same
idea.

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

## Current implication

This note records the preferred pattern for future API work:

- do not reintroduce wrapper-only execution layers
- keep construction on the structural noun rather than on the module
- keep execution methods type-level unless they own real runtime state
- keep examples pointed at the direct `State` / `Method` / `Evolution` seams

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
