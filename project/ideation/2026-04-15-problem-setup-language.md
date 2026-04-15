# Ideation: Problem Setup Language

Date: 2026-04-15

## Purpose

This note is a design artifact, not product documentation. The goal is to write
the *ideal* problem-setup surface for flux before further example cleanup,
benchmark expansion, or API growth. The sketches below are intentionally a bit
ahead of the current implementation.

The pressure test for each sketch is:

- Can a domain expert read the setup and recognize the physics?
- Are the structural nouns clear?
- Are cache and ownership boundaries explicit?
- Do the hot paths worth benchmarking follow naturally from the setup?

## Proposed top-level language

The setup surface should eventually make four structural objects explicit:

1. `Topology`
Represents incidence, orientation, and entity counts. This is where `d` lives.

2. `Geometry`
Represents embedding-dependent geometric data and mesh motion. For fixed meshes,
this may be little more than vertex positions plus derived measures. For moving
meshes, this becomes state.

3. `Metric`
Represents the inner-product structure used by metric-dependent operators such
as `★`, codifferentials, and Laplacians. Euclidean is the default structural
special case, not a boolean mode.

4. `Problem` or `System`
Owns the assembled operators, boundary conditions, and evolution/solve logic
for one PDE setup.

The current `Mesh` likely still contains both topology and some geometry. That
is acceptable today, but the setup language should stop pretending those are the
same thing.

## Sketch 1: User-defined fixed-mesh hyperbolic system

This is the baseline setup sketch. It should be the easiest problem to read and
the least magical. The important point is that flux should provide the language
for defining a system cleanly, not a built-in catalog of named PDE systems.

```zig
const flux = @import("flux");

const Topology = flux.topology.Mesh(2, 2);
var mesh = try Topology.rectangular(allocator, .{
    .cells_x = 64,
    .cells_y = 64,
    .width = 1.0,
    .height = 1.0,
});
defer mesh.deinit(allocator);

const geometry = flux.geometry.embedded(&mesh);
const metric = flux.geometry.Metric(.euclidean);

const MaxwellLikeSystem = flux.systems.DecSystem(Topology, struct {
    pub const State = struct {
        electric: flux.forms.Cochain(Topology, 1, flux.forms.Primal),
        magnetic: flux.forms.Cochain(Topology, 2, flux.forms.Primal),
    };

    pub fn assemble(
        allocator: std.mem.Allocator,
        mesh_ref: *const Topology,
        geometry_ref: anytype,
        metric_ref: anytype,
    ) !type {
        _ = allocator;
        _ = mesh_ref;
        _ = geometry_ref;
        _ = metric_ref;
        return struct {};
    }

    pub fn initialState(allocator: std.mem.Allocator, mesh_ref: *const Topology) !State {
        _ = allocator;
        _ = mesh_ref;
        return undefined;
    }

    pub fn rhs(self: *@This(), state: *const State) !State {
        _ = self;
        _ = state;
        return undefined;
    }
});

var system = try MaxwellLikeSystem.init(allocator, .{
    .mesh = &mesh,
    .geometry = &geometry,
    .metric = &metric,
    .boundary = .{ .essential = .all },
});
defer system.deinit();

var state = try system.initialState();
defer state.deinit();

var evolution = try flux.integrators.leapfrog.Evolution(@TypeOf(system)).init(allocator, .{
    .system = &system,
    .time_step = 1.0e-3,
    .observers = .{
        system.observe.scalar("energy"),
        system.observe.scalar("constraint"),
    },
});
defer evolution.deinit();

try evolution.advance(&state, .{ .steps = 1000 });
```

### What this sketch is trying to say

- `Mesh` is still the concrete mesh constructor, but the user sees separate
  `geometry` and `metric` objects immediately.
- The user chooses the mathematical family at the system level, but the PDE
  itself is user-defined rather than blessed by a built-in catalog of named
  equations.
- Boundary conditions belong to the problem setup, not to operator assembly in
  the hot path.
- The integrator consumes a system and a state. It should not need to know how
  the system assembled its operators.

### Nouns, verbs, adjectives

- Nouns: `Mesh`, `Geometry`, `Metric`, `System`, `State`, `Evolution`
- Verbs: `init`, `initialState`, `observe`, `advance`
- Adjectives: `Dec`, `euclidean`, `essential`

### Cache boundary

For this problem, the cache boundary should be:

- topology-owned: incidence
- geometry-owned: measures derived from embedding
- system-owned: assembled operators and boundary-conditioned variants
- evolution-owned: step counters, observer staging, transient integration state

### Benchmark targets implied by the setup

- repeated `d` application in the hyperbolic update
- repeated metric-dependent operator application if FEEC mass or `★` is used
- observer overhead relative to timestep cost
- one-time setup cost for assembling the system

## Sketch 2: User-defined metric-aware elliptic or parabolic system

This sketch exists to force metric semantics to be explicit and to test whether
the language scales from hyperbolic stepping to elliptic/parabolic solves.

```zig
const flux = @import("flux");

const Topology = flux.topology.Mesh(2, 2);

var mesh = try Topology.surfaceSphere(allocator, .{
    .refinement = 4,
});
defer mesh.deinit(allocator);

const geometry = flux.geometry.embedded(&mesh);
const metric = flux.geometry.Metric(.riemannian, .{
    .from_embedding = true,
});

const DiffusionLikeSystem = flux.systems.FeecSystem(Topology, struct {
    pub const State = flux.forms.feec.Form(
        flux.forms.feec.WhitneySpace(Topology, 0),
    );
});

var problem = try DiffusionLikeSystem.init(allocator, .{
    .mesh = &mesh,
    .geometry = &geometry,
    .metric = &metric,
    .boundary = .none,
});
defer problem.deinit();

var field = try problem.initialField(.{
    .analytic = diffusion_bump,
});
defer field.deinit();

var solver = try flux.solvers.implicit.backwardEuler(@TypeOf(problem)).init(allocator, .{
    .system = &problem,
    .time_step = 5.0e-4,
    .linear_solver = .{
        .cg = .{
            .preconditioner = .diagonal,
            .tolerance = 1.0e-10,
        },
    },
});
defer solver.deinit();

try solver.advance(&field, .{ .steps = 200 });
```

### What this sketch is trying to say

- Metric choice is visible at setup time.
- FEEC is selected because the weak-form and mass/stiffness structure matter.
- The solver stack should read as part of the PDE language, not as detached
  sparse-matrix plumbing.
- If the metric changes, users should know which operators are invalidated.

### Nouns, verbs, adjectives

- Nouns: `Mesh`, `Geometry`, `Metric`, `System`, `Field`, `Solver`
- Verbs: `init`, `initialField`, `advance`
- Adjectives: `Feec`, `riemannian`, `diagonal`

### Cache boundary

For this problem, the cache boundary should make invalidation obvious:

- changing topology invalidates everything
- changing geometry may invalidate embedding-derived measures
- changing metric invalidates metric-dependent FEEC operators
- changing boundary conditions invalidates only the boundary-conditioned solve
  structure, not topology itself

### Benchmark targets implied by the setup

- assembly of mass and stiffness operators
- repeated sparse solve cost per timestep
- preconditioner effectiveness on user-visible wall time
- metric-update cost if geometry or metric is changed between solves

## Sketch 3: User-defined moving-geometry or ALE-flavored transport

This sketch is intentionally horizon-facing. Its job is to force the language to
separate topology from evolving geometry before the implementation needs it.

```zig
const flux = @import("flux");

const Topology = flux.topology.Mesh(2, 2);

var mesh = try Topology.rectangular(allocator, .{
    .cells_x = 96,
    .cells_y = 48,
    .width = 2.0,
    .height = 1.0,
});
defer mesh.deinit(allocator);

var geometry = try flux.geometry.MovingGeometry(Topology).init(allocator, .{
    .mesh = &mesh,
    .embedding = .current_positions,
});
defer geometry.deinit();

var metric = flux.geometry.Metric(.euclidean);

const AleTransportSystem = flux.systems.AleSystem(Topology, struct {
    pub const State = struct {};
});

var problem = try AleTransportSystem.init(allocator, .{
    .mesh = &mesh,
    .geometry = &geometry,
    .metric = &metric,
    .mesh_velocity = prescribed_mesh_velocity,
});
defer problem.deinit();

var state = try problem.initialState(.{
    .density = gaussian_blob,
    .velocity = rigid_rotation,
});
defer state.deinit();

try problem.advanceGeometry(.{
    .time_step = 1.0e-3,
    .steps = 500,
    .state = &state,
});
```

### What this sketch is trying to say

- Geometry is now plainly stateful and evolves in time.
- The topology remains fixed unless explicitly remeshed.
- The system can ask for reassembly of only the metric-dependent pieces after a
  geometry update.
- This should be possible without turning the mesh itself into a mutable bag of
  cached operators.

### Nouns, verbs, adjectives

- Nouns: `Mesh`, `MovingGeometry`, `Metric`, `System`, `State`
- Verbs: `init`, `initialState`, `advanceGeometry`
- Adjectives: `Ale`, `moving`

### Cache boundary

This is the most important sketch for future-proofing:

- topology-only operators should survive geometry motion untouched
- geometry-dependent derived quantities should be recomputed incrementally if
  possible
- metric-dependent assembled operators should be invalidated explicitly and only
  when needed

### Benchmark targets implied by the setup

- geometry-update cost
- operator invalidation and reassembly cost after geometry motion
- ratio of pure topology work to metric/geometry recomputation
- total cost per coupled PDE step

## Design pressure revealed by the sketches

### 1. The family choice should likely sit at the problem or system layer

Users should not have to manually assemble a bag of DEC and FEEC operators just
to express a standard problem. The problem setup should name the mathematical
family directly, while still allowing the actual PDE system to be defined by
the user rather than by a built-in framework catalog.

### 2. `Metric` should be explicit, but probably not alone

The sketches suggest `Metric` is necessary but not sufficient. A broader
`Geometry` noun likely needs to exist so the API can distinguish:

- fixed embedding-derived geometry
- prescribed metric on fixed topology
- time-evolving geometry with possibly evolving metric dependence

### 3. Caches should remain outside topology ownership

All three sketches argue against moving assembled operators back onto the mesh.
What users want is not mesh-owned numerics, but a setup surface that *reads*
coherently while keeping ownership boundaries honest.

### 4. Benchmark design should derive from setup and execution language

The benchmark suite should answer user-relevant questions:

- what costs dominate setup?
- what costs dominate one step?
- what costs dominate repeated solve/apply loops?
- which costs change when geometry or metric changes?

This is stronger than “benchmark all public operators.” It connects the suite to
real usage rather than implementation accidents.

## Recommendations

1. Treat this note as input to example redesign, not as documentation to make
   current examples mimic mechanically.

2. Use these sketches to draft a narrower design note for the stable top-level
   nouns:
   `Mesh`, `Geometry`, `Metric`, `System`, `State`, `Evolution` or `Solver`.

3. When designing the abstraction-audit skill, require it to check whether new
   APIs move the project closer to these nouns by deleting weak layers rather
   than wrapping them.

4. Derive future benchmark additions from the hot paths named in these sketches.

## Open questions

- Should `Geometry` be a first-class public noun now, or remain implicit until a
  moving-geometry use case lands?
- Where should mixed FEEC/DEC formulations appear in the public setup surface:
  separate `MixedSystem`, explicit expert composition, or a problem-specific
  constructor?
- Should integrator-facing types be called `Evolution`, `Stepper`, `Solver`, or
  something problem-class-specific?
- Which parts of boundary condition handling belong to `System` versus lower
  operator layers?
