# Architecture Execution

Date: 2026-04-15

## Purpose

This note defines the durable execution language for flux: what owns repeated
time stepping, what contract a stepper must satisfy, and which wrapper patterns
should be deleted rather than normalized.

This is the canonical reference for execution-side API reviews. If a proposed
execution type is hard to justify against this note, it should probably not
exist.

## Core rule

Execution should read in terms of a small set of structural roles:

- `System` owns assembled operators, solve state, boundary-conditioned variants,
  and other problem-local caches
- `State` is the evolving unknown
- `Evolution` owns repeated-step orchestration
- `TimeStepper` is a comptime contract, not a forwarding wrapper

If examples need builder nouns, adapter nouns, or wrapper-only tests to fit the
execution layer, the seam is wrong.

## Preferred execution boundary

The default shape is:

- define a comptime stepper concept
- validate that concept at the `Evolution` boundary
- pass the conforming stepper directly

The default shape is **not**:

- define a stepper concept
- wrap it in a `TimeStepper(...)` forwarding type
- wrap that again for a runner or evolution helper
- add builder objects whose only job is to manufacture the wrapper

In other words: the boundary is the concept check.

## `TimeStepper` is a concept

`TimeStepper` should mean the direct execution contract required by
`Evolution`. It should not also name a second wrapper noun.

The preferred contract is intentionally small:

- `pub fn step(self: *Stepper, allocator: std.mem.Allocator) !void`
- `pub fn deinit(self: *Stepper, allocator: std.mem.Allocator) void`

If execution later needs stronger guarantees, extend the concept. Do not add a
wrapper type whose only job is to re-check the same declarations and forward
them.

This implies:

- delete `src/time_stepper.zig` as a wrapper home
- keep the concept near the execution layer, ideally as part of the evolution
  module surface
- make `Evolution` validate `TimeStepper` directly

## `Evolution` is the execution noun

`Evolution` earns its existence when it owns real orchestration state such as:

- allocator and lifetime management for execution-owned state
- step counts
- run timing
- listener dispatch
- exact/reference auxiliary state
- snapshot or observer staging

That is a real noun. It should stay.

What should **not** exist between `Evolution` and the stepper:

- forwarding wrappers
- builder objects with a single `initStepper` method
- adapter types that only rename `step` and `deinit`

If setup-time seeding is required before stepping begins, that belongs in one
of two places:

- system initialization, if it is inherent to system construction
- stepper initialization, if it is specific to one execution story

It should not force a second builder noun by default.

## Relationship to `System`

`System` and `Evolution` should remain separate.

- `System` is the PDE/model noun
- `Evolution` is the repeated-step execution noun

The right simplification is a more direct contract between them, not merging
them and not inserting additional wrapper layers.

A stepper may hold:

- a pointer to system-owned state
- references to state storage
- execution-local counters or scratch

But those choices should be expressed directly in the stepper type the caller
passes to `Evolution`.

## Example pressure tests

Good execution design should make these examples small:

1. an explicit hyperbolic system with a custom stepper
2. an implicit parabolic system whose stepper seeds a linear solve state once
3. a reference or exact-solution run with observer/listener hooks

If one of these requires:

- `HeatStepperBuilder`
- `SurfaceStepperBuilder`
- `StepperWrapper`
- `TimeStepper(Strategy)`

then the execution seam is too indirect.

## Test policy

Prefer tests that prove:

- `Evolution` rejects a non-conforming stepper at comptime
- a conforming stepper advances state correctly
- listeners and exact/reference hooks observe the right run events
- PDE invariants hold over repeated steps

Do not preserve tests whose only purpose is:

- “the wrapper forwards to the inner type”
- “the builder returns the wrapper”

Those tests protect ceremony, not behavior.

## Deletion guidance

The default cleanup order for execution abstraction drift is:

1. delete forwarding wrappers
2. inline one-shot builder setup into direct init
3. make `Evolution` consume the direct concept
4. only then ask whether any additional noun is still necessary

If a proposed fix adds execution code, the burden is to explain why deletion
and direct concept validation are insufficient.

## Current intended direction

The current intended direction is:

- `TimeStepper` becomes a comptime concept only
- `Evolution` remains the execution owner
- `Evolution` validates the direct stepper contract
- `src/time_stepper.zig` should disappear as a wrapper-focused module
- example-local stepper builders and adapter wrappers should be removed as the
  execution seam is tightened

This note records the target shape. Individual code changes should converge
toward it incrementally, but they should not add new wrapper layers in the
meantime.
