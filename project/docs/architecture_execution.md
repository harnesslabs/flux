# Execution Architecture

Date: 2026-04-17

## Purpose

This note defines the execution-side nouns for flux and the preferred
construction pattern for time evolution.

The goal is a small, explicit execution language:

- `System` owns one runtime PDE instance
- `Method` names a time-integration family
- `Evolution` owns repeated stepping, time, counters, and configured listeners
- `Measurement` names optional observable quantities exposed by a system
- reference/comparison helpers stay outside core execution

## Core nouns

### `System`

`System` is the owned runtime PDE instance for one problem family.

It may own:

- cochains and field storage
- operator contexts and assembled operators
- mesh references
- boundary-condition data
- assembled linear systems
- problem-specific caches
- counters intrinsic to the physical update

It should not be wrapped just to satisfy the execution API.

Systems expose:

- structural fields through `pub const Field`
- optional scalar observables through `pub const Measurement`
- optional integration contracts that method families validate at comptime

### `Method`

`Method` is a time-integration family: leapfrog, forward Euler, backward
Euler, or a future adaptive controller.

Default rule:

- keep methods as comptime families or comptime-specialized types
- do not invent a runtime stepper object that re-stores `system`, `dt`, or
  other execution-owned state

Method families may expose:

- `pub fn advance(allocator, system, dt) !void`
- `pub const Options = struct { ... }`
- `pub fn initialize(allocator, system, dt, options) !void`

Method-specific configuration flows through `Evolution.Config.methodOptions(...)`.

### `Evolution`

`Evolution` is the execution owner.

It owns:

- the primary system handle passed to it
- `dt`
- configured run length (`steps`)
- current time
- step count
- configured listeners
- stored method options

It does **not** own a separate runtime stepper object by default.

### `Measurement`

`Measurement` is a system-declared observable quantity that is not a structural
field, such as total energy.

Execution listeners may request fields and measurements from a system.
Comparison- or refinement-driven quantities may later be expressed through the
same measurement surface once that API settles.

## Construction pattern

The standard construction path is now:

```zig
var evolution = try flux.evolution.Evolution(SystemType, flux.integrators.Leapfrog)
    .config()
    .dt(dt)
    .steps(steps)
    .listen(flux.listeners.Progress(writer))
    .init(allocator, system);
```

Interpretation:

- `Evolution(System, MethodFamily)` is the type-level execution noun
- `.config()` begins the staged configuration phase
- `.dt(...)` configures execution-owned timestep policy
- `.steps(...)` configures run length
- `.listen(...)` attaches event listeners that live on the evolution object
- `.methodOptions(...)` configures method-specific settings when needed
- `.init(...)` constructs the runtime execution object
- `.run()` executes the configured run using the stored listeners

Method families are specialized internally by `Evolution`; users should not
spell `Leapfrog(System)` separately at ordinary call sites.

## Ownership rules

### `dt`

`dt` belongs to `Evolution`, not to a thin method wrapper.

Reason:

- time-step policy is part of execution
- fixed-step and future adaptive-step methods should share one construction
  surface
- method-specific controller settings can stay downstream in
  `.methodOptions(...)`

### Listeners

Listeners belong to `Evolution`, not to the shared example runner.

Reason:

- they are execution attachments
- they observe run begin/step/end events
- they should be configured once and reused by `run()`

### Measurements

Measurements belong to the `System` surface.

Reason:

- canonical observables such as total energy are properties of the chosen PDE
  system, not of execution itself
- listeners can stay generic if systems expose a small measurement interface

## Wrapper admission rule

Do not introduce a runtime `Stepper`/`Integrator` object unless it owns a real
runtime identity such as:

- mutable controller state that should not live on `Evolution`
- owned scratch or cached factorization state that should not live on `System`
- a representation boundary that materially simplifies the public language

If the proposed object only:

- holds `state`
- holds `system`
- holds `dt`
- forwards one `step` call

then it should not exist.

## Current implication

The default flux execution pattern is now:

- `Evolution(System, MethodFamily).config().dt(...).steps(...).listen(...).init(...)`
- methods are type-level by default
- `Evolution` specializes method families internally
- method options flow through `Evolution.Config`
- listeners attach during configuration and run through `evolution.run()`
- systems expose structural fields and optional measurements
- reference/comparison helpers remain outside core execution for now

This note should be updated before introducing adaptive execution features,
public comparison/probe APIs, or any new public execution noun.
