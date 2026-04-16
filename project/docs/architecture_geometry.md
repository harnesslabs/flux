# Architecture Geometry

Date: 2026-04-15

## Purpose

This note defines the intended boundary between topology, geometry, and metric
structure in flux.

The implementation is still earlier than the final public language in some
places. That is acceptable. What is not acceptable is letting temporary module
placement harden into the long-term abstraction language by accident.

## Core rule

Users should be able to distinguish:

- what is purely topological
- what depends on embedding geometry
- what depends on metric structure
- what is system-owned assembly derived from those ingredients

If those roles blur together at the API boundary, later generalization to
Riemannian metrics, moving geometry, or alternate assembly families becomes
needlessly expensive.

## Preferred nouns

### `Mesh`

`Mesh` is the concrete discrete domain object used today.

At present it still carries:

- topology
- vertex positions
- some embedding-derived geometric data
- some precomputed operator-support data

That mixed reality is tolerated as an implementation stage, not celebrated as
the final public language.

### `Geometry`

`Geometry` is the noun for embedding-dependent structure:

- coordinates
- edge lengths
- face areas
- cell volumes
- future mesh motion or geometry updates

Even when there is not yet a separate runtime `Geometry` object at every call
site, this noun should guide design decisions. Geometry is not merely “whatever
is stored on mesh today.”

### `Metric`

`Metric` is the noun for inner-product structure:

- Euclidean as the explicit default specialization
- Riemannian as a genuine extension

The important point is conceptual: metric choice is a structural input, not an
operator-local policy flag.

That means the long-term language should read like:

- `Metric(.flat)` or `Metric(.euclidean)`
- `Metric(.riemannian, ...)`

and not like a hidden mode buried inside one operator family.

## Ownership boundary

The intended ownership split is:

- topology-owned: incidence, orientations, entity counts
- geometry-owned: embedding-derived measures and future geometric updates
- metric-owned: inner-product data and metric-specialized structure
- system-owned: assembled operators, boundary-conditioned variants, and caches

This is why mesh-owned operator caches are a bad default public model even when
some low-level data is currently stored on `Mesh`.

## Current implementation status

Today, metric-aware execution is only partially surfaced:

- metric-aware Hodge star support currently lives in
  `src/operators/hodge_star.zig`
- the flat specialization remains the default operator-context path
- there is not yet a first-class top-level geometry namespace that owns
  `Metric`

This is a temporary implementation shape, not the desired final language.

## Rule for API growth during the transition

Until geometry and metric become more explicit public modules:

- do not add new public APIs that make metric look like an implementation
  detail of Hodge star alone
- do not add new “mode” booleans when the real noun is metric structure
- do not treat mesh storage location as proof of long-term ownership

New code may temporarily use the existing operator-local metric types, but new
design should speak in terms of the geometry/metric boundary above.

## Family interaction

The topology/geometry/metric split is orthogonal to DEC versus FEEC family
choice.

- `d` is topology-driven
- Hodge stars, codifferentials, Laplacians, and weak inner products are metric
  dependent
- system-level family choice should decide which assembled operators are used,
  not collapse geometry and metric into one opaque implementation switch

Family adjectives refine operator and system assembly. They do not replace the
underlying structural nouns.

## Invalidation guidance

The intended invalidation logic is:

- topology changes invalidate everything
- geometry changes invalidate embedding-derived measures and dependent assembly
- metric changes invalidate metric-dependent operators and solves
- boundary-condition changes invalidate only the affected system-owned
  conditioned structures

This should remain legible at the API level. If a caller cannot tell what a
metric change invalidates, the abstraction boundary is too weak.

## Pressure tests

A good geometry boundary should support:

1. fixed Euclidean meshes
2. fixed meshes with user-supplied Riemannian metric data
3. future moving-geometry or ALE-flavored problems

If any of these would force a new wrapper family or another operator-local
special case API, the geometry language is still too accidental.

## Current intended direction

The current intended direction is:

- keep flat-specialized operator contexts working as the zero-cost default
- treat operator-local `Metric` support as a temporary staging point
- move toward a public language where geometry and metric are structural setup
  inputs rather than hidden operator options

This note does not force an immediate refactor of all operator code. It does
set the standard for future API decisions so temporary placement does not
become permanent design by inertia.
