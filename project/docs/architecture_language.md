# Architecture Language

Date: 2026-04-15

## Purpose

This note defines the preferred public language for setting up and executing
problems in flux. It is the canonical destination for the durable parts of the
problem-setup ideation work.

The goal is not to freeze every API detail today. The goal is to make the
structural nouns, verbs, and qualifiers explicit so new APIs converge toward one
coherent language instead of growing ad hoc.

## Core rule

Users should be able to read a flux setup and recognize:

- what is topology
- what is geometry
- what is metric structure
- what is the PDE system
- what is the execution or solve mechanism

If a setup reads like a pile of helper types, wrappers, builders, or
case-specific names, the language is still wrong.

## Preferred top-level nouns

### `Mesh`

`Mesh` is the concrete discrete domain constructor and storage object.

Today it may still carry both topology and some geometry-derived data. That is
acceptable as an implementation reality, but the public language should not
pretend those concerns are fundamentally identical.

`Mesh` should stay a structural noun, refined by adjectives or comptime
parameters such as:

- embedding dimension
- topological dimension
- family of canonical geometry constructor

### `Geometry`

`Geometry` is the noun for embedding-dependent geometric structure and, later,
mesh motion.

This layer should eventually make it possible to distinguish:

- fixed embedding-derived geometry
- prescribed geometry on fixed topology
- moving geometry whose state evolves in time

The existence of this noun matters even before a fully separate runtime object
exists everywhere in the implementation. It is the right mental boundary for
future API design.

### `Metric`

`Metric` is the noun for the inner-product structure used by metric-dependent
operators such as `★`, codifferentials, and Laplacians.

The default special case should be explicit:

- `Metric(.euclidean)`
- `Metric(.riemannian, g)`

This is preferred over vague “mode” flags because it names the mathematical
object rather than the policy choice.

`Metric` is likely not the whole geometry story by itself. It should be treated
as one structural layer inside a broader geometry boundary.

### `System`

`System` is the noun for one PDE setup over a mesh, geometry, and metric.

The system owns:

- the chosen mathematical family or families
- assembled operators and boundary-conditioned variants
- PDE-specific state shape
- problem-specific setup semantics

Flux should provide the language for defining systems cleanly. It should not
grow into a built-in catalog of named PDE APIs such as “Maxwell system,”
“Euler system,” and so on as first-class core-library nouns.

The library should instead let users define those systems against a stable
structural language.

### `State`

`State` is the evolving or solved-for discrete unknown associated with a
`System`.

It should be named by mathematical role, not by framework plumbing. For
example, a state may contain fields, forms, or coefficients, but it should not
exist only to satisfy a wrapper seam.

### `Evolution` / `Solver`

Execution nouns should reflect the actual job being done:

- `Evolution` for repeated time stepping
- `Solver` for linear, nonlinear, or implicit solve processes

Do not invent extra nouns when a strong verb on an existing structural noun
would do. If an execution object exists, it should own real orchestration state,
not only forward to another type.

`Evolution` should usually be a capability over a `System`, not the same noun.
The separation matters:

- `System` defines the PDE/model structure and state semantics
- `Evolution` defines one repeated-step execution story over systems that
  support stepping

Some systems do not evolve in time at all. Others may admit multiple execution
stories (explicit, implicit, split, multirate) over the same underlying system.
For that reason, the project should resist collapsing `System` and `Evolution`
into one default public noun.

The tighter direction is:

- keep `System` and `Evolution` separate
- make `Evolution` depend on a direct evolvable-system contract
- avoid extra wrapper nouns between the system, the stepper concept, and the
  orchestration object

## Family choice

The mathematical family choice should sit at the system layer.

Users should not have to manually assemble a bag of DEC and FEEC operators just
to express a standard problem. The setup should name the family directly, while
still allowing:

- user-defined PDE systems
- mixed FEEC/DEC formulations when the semantics genuinely differ
- explicit expert composition when a problem really needs it

The project should resist carrying parallel public APIs that express the same
operation twice merely because one path is older or more convenient.

## Preferred verbs

The language should lean on a small set of strong verbs:

- `init`
- `assemble`
- `advance`
- `solve`
- `observe`
- `project`
- `lift`

These verbs should sit on stable nouns. Avoid creating a new noun when the real
operation is simply “do X to this object.”

For FEEC weak operators specifically, prefer one shared `assemble` seam for
local-to-global scatter, with family-specific local kernels supplying the
simplex-local mathematics. Do not duplicate orientation/scatter plumbing inside
each weak operator family helper.

## Qualifiers should usually be adjectives

Many distinctions in flux are qualifiers, not reasons to mint fresh public
nouns. Prefer adjectives or comptime parameters for:

- dimension
- degree
- primal versus dual
- scalar type
- DEC versus FEEC family choice
- Euclidean versus Riemannian metric choice

Create a new noun only when the underlying mathematical or computational role is
different, not merely because one more policy axis was introduced.

## Ownership and cache boundaries

The language should make the following ownership boundaries legible:

- topology-owned: incidence and purely topological structure
- geometry-owned: embedding-derived measures and geometric updates
- metric-owned: metric-dependent structure
- system-owned: assembled operators, boundary-conditioned variants, and
  problem-local caches
- execution-owned: repeated-step orchestration, solver state, observer staging

This is the main reason the project should avoid mesh-owned operator caches as a
default public model. What users want is a coherent setup language, not blurred
ownership.

## Example pressure tests

A good setup language should support at least these three shapes cleanly:

1. User-defined fixed-mesh hyperbolic system
2. User-defined metric-aware elliptic or parabolic system
3. User-defined moving-geometry or ALE-flavored system

If one of those forces a pile of builders, wrappers, or example-only helper
nouns, the language needs revision.

## Benchmark implication

Benchmarks should be chosen from the hot paths implied by this language:

- setup and assembly cost
- repeated operator application cost
- repeated solve cost
- observer overhead
- invalidation and reassembly cost under geometry or metric changes

This is stronger than benchmarking whichever functions happen to exist. The
benchmark suite should measure user-relevant execution paths.

## API review checklist

Before adding a new public type or helper, ask:

1. Is this a stable structural noun?
2. Could this instead be an adjective on an existing noun?
3. Could this instead be a verb on an existing noun?
4. Does this clarify topology, geometry, metric, system, state, or execution?
5. Would the new shape make example setups smaller and clearer?

If the answer is weak, do not add the new noun.

## Open edges

The following points are intentionally not frozen yet:

- whether `Geometry` should become a first-class runtime object immediately
- where mixed FEEC/DEC formulations should appear in the public setup surface
- the final naming split between `Evolution`, `Stepper`, and `Solver`
- the exact direct contract between `System` and `Evolution`
- the exact boundary between system-level and lower-level boundary-condition APIs

These should evolve under this language, not around it.
