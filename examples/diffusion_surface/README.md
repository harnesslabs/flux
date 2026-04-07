# Diffusion On A Curved Surface

This example solves

\[
\partial_t u = \Delta_{S^2} u
\]

on a triangulated unit sphere. The spatial operator is assembled from the
Riemannian DEC ingredients on a 2D reference mesh:

\[
S = d_0^\top \star_1^{(g)} d_0,
\qquad
(M + \Delta t\,S)u^{n+1} = M u^n.
\]

The exact mode used for verification is the degree-1 spherical harmonic
`u(x, y, z, t) = e^{-2t} z`, since `Δ_{S^2} z = 2 z`.

## Run

```sh
zig build -Doptimize=ReleaseFast run-diffusion-surface -- --refinement 3 --frames 8
zig build -Doptimize=ReleaseFast run-diffusion-surface -- --help
```

## Output

Snapshots are written to `output/diffusion_surface/` by default:

- `diffusion_surface_XXXX.vtu` stores point fields `temperature`, `temperature_exact`, and `temperature_error`
- `diffusion_surface.pvd` indexes the time series for ParaView

## Verification

The example test suite verifies:

- spherical-harmonic decay error decreases on refinements `1`, `2`, and `3`
- the finest tested refinement stays below `1e-3` in normalized lumped `L²`
- the standalone binary configuration runs through the public package entrypoint
