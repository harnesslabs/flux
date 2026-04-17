# Component Scope Map

This file defines the codebase's component boundaries. Skills use it to limit
context window usage — when working on an issue in a given domain, load only
that component's files and its direct dependencies.

## Components

### topology
**Domain label:** `domain/topology`
**Owns:** `src/topology/`
**Dependencies:** `src/math/` (sparse matrix types)
**Description:** CW complex representation, incidence/boundary matrices, mesh
construction, dual mesh computation, geometric data (edge lengths, face areas,
dual volumes).

### forms
**Domain label:** `domain/forms`
**Owns:** `src/forms/`
**Dependencies:** `src/topology/` (mesh types for comptime parameterization)
**Description:** Discrete k-forms (cochains) parameterized by mesh and degree at
comptime. Arithmetic operations on cochains (add, scale, negate).

### operators
**Domain label:** `domain/operators`
**Owns:** `src/operators/`
**Dependencies:** `src/topology/`, `src/forms/`, `src/math/`
**Description:** Discrete exterior derivative, Hodge star, Laplacian, operator
composition. Each operator is an independent, testable unit.

### concepts
**Domain label:** `domain/topology`
**Owns:** `src/concepts/`
**Dependencies:** `src/topology/`
**Description:** Comptime concept validators (`MeshConcept`). Enforce compile-time
contracts on user-defined types so alternative mesh implementations (mesh views,
imported meshes, n-dimensional meshes) are guaranteed compatible with the operator stack.

### integrators
**Domain label:** `domain/operators` (shares label with spatial operators)
**Owns:** `src/integrators/`
**Dependencies:** none (generic — parameterized on system types at comptime)
**Description:** Execution and time-integration patterns: `Evolution` owns
runtime time management and listeners, while integrator methods such as
leapfrog and forward Euler are parameterized at comptime on system types.
Reference studies are separate higher-order helpers rather than built into
`Evolution`. Concrete physics modules provide the update functions; the
execution layer configures and applies them without wrapper-only runtime
steppers.

### math
**Domain label:** `domain/build` (no dedicated label yet)
**Owns:** `src/math/`
**Dependencies:** none
**Description:** Low-level math primitives — sparse matrix types (CSR), linear
algebra utilities. Foundation layer with no upward dependencies.

### io
**Domain label:** `domain/io`
**Owns:** `src/io/`
**Dependencies:** `src/topology/`, `src/forms/`
**Description:** VTK export for visualization. Reads mesh and cochain data,
writes `.vtu` files.

### examples
**Domain label:** `domain/em`
**Owns:** `examples/maxwell/`, `examples/euler/`, `examples/diffusion/`, `examples/common/`, `examples/commands.zig`, `examples/app.zig`, `examples/acceptance.zig`
**Dependencies:** flux library (via package import)
**Description:** End-to-end physics examples and their shared CLI/output
infrastructure. Maxwell and Euler are dimension-dispatched from one module per
physics family; diffusion dispatches plane versus sphere from one module.

### cli
**Domain label:** `domain/build`
**Owns:** `src/main.zig`
**Dependencies:** none (minimal stub)
**Description:** Placeholder CLI entry point. Directs users to examples.

### library
**Domain label:** `domain/build`
**Owns:** `src/root.zig`
**Dependencies:** all `src/` modules
**Description:** Library entry point — re-exports public API. Changes here
affect what downstream consumers see.

### ci
**Domain label:** `domain/ci`
**Owns:** `.github/`, `build.zig`
**Dependencies:** none (build-system level)
**Description:** CI workflows, GitHub settings, build configuration, label
definitions, issue templates.

### project
**Domain label:** `domain/docs`
**Owns:** `project/`, `.claude/`
**Dependencies:** none
**Description:** Epoch planning, decision logs, retrospectives, skill
definitions, vision and horizons documents.

## Dependency graph

```mermaid
flowchart TD
    math --> topology
    topology --> forms
    math --> operators
    topology --> operators
    forms --> operators
    topology --> io
    forms --> io

    subgraph "flux library (src/)"
        math
        topology
        forms
        operators
        io
        integrators
        concepts
    end

    subgraph "examples/"
        examples_suite["maxwell / euler / diffusion"]
    end

    flux_library["flux (package)"] --> examples_suite
```

## Usage

When `/tackle` identifies an issue's domain, it should:
1. Look up the component in this file
2. Read the **Owns** files and **Dependencies** files
3. Do not read files outside this set unless the issue explicitly crosses component boundaries
4. If a cross-component change is needed, state which components are involved and why
