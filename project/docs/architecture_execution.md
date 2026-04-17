# Execution Architecture

Date: 2026-04-17

## Purpose

This note defines the execution-side nouns for flux and the preferred
construction pattern for time evolution.

The goal is a small, explicit execution language:

- `State` owns the evolving PDE data and problem-owned caches
- `Method` names the time-integration algorithm
- `Evolution` owns repeated stepping, time, counters, listeners, and optional
  exact/reference state

## Core nouns

### `State`

`State` is the simulation state for one problem family.

It may own:

- cochains and field storage
- operator contexts and assembled operators
- mesh references
- problem-specific caches
- counters that are intrinsic to the physical update itself

It should not be wrapped just to satisfy the execution API.

### `Method`

`Method` is the integration algorithm: leapfrog, forward Euler, backward Euler,
 or a future adaptive method.

Default rule:

- if the method has no real runtime identity, keep it as a comptime type with
  `advance(...)`
- do not invent a runtime `Stepper` object that re-stores `state`, `dt`, or
  other execution-owned data

Methods may expose:

- `pub fn advance(allocator, state, dt) !void`
- `pub const Options = struct { ... }`
- `pub fn initialize(allocator, state, dt, options) !void` when construction
  needs method-specific setup

Method options are configuration, not a second execution noun.

### `Evolution`

`Evolution` is the execution owner.

It owns:

- the primary state reference/value passed to it
- `dt`
- current time
- step count
- listeners and run orchestration
- optional exact/reference auxiliary state
- stored method options

It does **not** own a separate runtime stepper object by default.

## Construction pattern

The standard construction path is:

```zig
const Evolution = flux.evolution.Evolution(StateType, MethodType, AuxType);

var evolution = try Evolution
    .config()
    .dt(dt)
    .methodOptions(.{ .source = source })
    .init(allocator, state, aux);
```

Interpretation:

- `Evolution(...)` is the type-level noun, like `Mesh(...)`
- `.config()` begins a staged configuration phase
- `.dt(...)` configures execution-owned timestep policy
- `.methodOptions(...)` configures method-specific settings when needed
- `.init(...)` constructs the runtime execution object

Examples should prefer the inline form above over:

```zig
const EvolutionType = flux.evolution.Evolution(...);
var evolution = try EvolutionType.config()...init(...);
```

unless the type name itself is needed later in type position.

## Ownership rules

### `dt`

`dt` belongs to `Evolution`, not to a thin method wrapper.

Reason:

- time-step policy is part of execution
- fixed-step and adaptive-step methods should share one construction surface
- future adaptive methods can refine how `Evolution` chooses or updates `dt`
  without changing the top-level noun split

### Method-specific setup

If a method needs setup during construction, that setup happens downstream of
`Evolution.Config.init(...)` through an optional method hook such as
`initialize(...)`.

Example:

- backward Euler methods may seed a linear-system solution buffer from the
  current state
- a future adaptive method may validate tolerances or initialize controller
  state

This keeps the construction surface flat while still allowing method-specific
 setup.

## Wrapper admission rule

Do not introduce a runtime `Stepper`/`Integrator` object unless it owns a real
runtime identity such as:

- mutable controller state that is separate from `Evolution`
- owned scratch or cached factorization state that should not live on `State`
- a representation boundary that materially simplifies the public language

If the proposed object only:

- holds `state`
- holds `dt`
- forwards one `step` call

then it should not exist.

## Example quality signal

Example code should read like:

- choose a `State`
- choose a `Method`
- configure `Evolution`
- run

When examples start inventing family-local wrapper structs just to satisfy the
execution interface, the interface is wrong.

## Current implication

The default flux execution pattern is now:

- `Evolution(State, Method, Aux).config().dt(...).init(...)`
- methods are type-level by default
- method options flow through `Evolution.Config`
- `Evolution` owns time management and run bookkeeping

This note should be updated before introducing adaptive execution features,
observer attachment changes, or any new public execution noun.
