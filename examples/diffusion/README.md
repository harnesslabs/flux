# Diffusion

Scalar diffusion examples with one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --help
zig build -Doptimize=ReleaseFast run-diffusion -- --surface sphere --help
```

## Surfaces

- `--surface plane` solves the implicit heat equation on the unit square with
  homogeneous Dirichlet boundary data.
- `--surface sphere` solves diffusion on a triangulated unit sphere using the
  degree-1 spherical harmonic `u(x, y, z, t) = e^{-2t} z` as the exact mode.

Both paths share the same top-level diffusion command and convergence-test
surface.

## Run

```sh
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --grid 32 --frames 8
zig build -Doptimize=ReleaseFast run-diffusion -- --surface sphere --refinement 3 --frames 8
```

## Verification

- The plane example verifies zero-data preservation, boundary pinning, and
  second-order spatial convergence.
- The sphere example verifies monotone error reduction under refinement and a
  tight finest-grid analytic error bound.
