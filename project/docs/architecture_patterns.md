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

- keep one direct stepper concept at the execution boundary
- let `Evolution` validate that concept directly
- pass the conforming stepper value directly
- expose an inferred value constructor for ordinary call sites

In other words, the interface boundary is the concept check, not a wrapper
around the concept.

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

## Pattern: keep type factories available, but do not teach them as the default constructor path

Many flux nouns are legitimately parameterized at comptime. That means a
type-level factory such as:

- `Mesh(2, 2)`
- `Cochain(MeshType, k, Duality)`
- `Evolution(StateType, StepperType, AuxType)`

may still be the right implementation tool.

But the existence of a type factory does **not** mean ordinary callers should
be forced to construct values through:

- `Foo(T, U, V).init(...)`

when the needed type parameters can be inferred from the value arguments.

The preferred public construction pattern is:

- keep the type factory for type-position use
- expose one value-level constructor or helper when callers mainly want a value
- teach the value-level path as the default usage

### Why

If a caller writes:

- `var evolution = flux.evolution.init(allocator, state, stepper, aux);`

the resulting value still has the full concrete type. We have not lost type
information. We have only stopped making the caller spell the same structure
twice.

This preserves the one-obvious-way rule better than teaching both:

- `const EvolutionType = flux.evolution.Evolution(State, Stepper, Aux);`
- `var evolution = EvolutionType.init(...)`

as equally normal construction patterns.

The explicit type factory remains useful in known type-position cases:

- return types
- field types
- comptime reflection
- nested declarations such as event payload types

Those are real use cases. They are not a reason to make the explicit type
constructor the default value-construction story.

### Anonymous literal caution

When a value-level constructor relies on type inference, callers should avoid
passing anonymous struct literals when a named runtime type is semantically
important.

Weak:

- `evolution.init(allocator, state, .{ .state = &state, .dt = dt }, aux)`

Stronger:

- `const stepper = EulerStepper{ .state = &state, .dt = dt };`
- `var evolution = evolution.init(allocator, state, stepper, aux);`

This keeps the call site explicit about the intended noun while still avoiding
manual `@TypeOf(...)` plumbing.

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
- do not teach type-factory-plus-`init` as the default value-construction path
  when a direct constructor can infer the type parameters
- keep examples pointed at the direct concept and value-construction seams

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
