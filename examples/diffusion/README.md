# Diffusion

One family module, one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --help
zig build -Doptimize=ReleaseFast run-diffusion -- --surface sphere --help
```

## Surfaces

| Selector | Purpose | Main knobs |
|---|---|---|
| `--surface plane` | implicit heat equation on the unit square | `--grid`, `--dt-scale`, `--steps` |
| `--surface sphere` | diffusion of the spherical harmonic `e^{-2t} z` | `--refinement`, `--final-time`, `--steps` |

## Run

```sh
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --grid 32 --frames 8
zig build -Doptimize=ReleaseFast run-diffusion -- --surface sphere --refinement 3 --frames 8
```

## Model

- Plane mode uses backward Euler on the square with homogeneous Dirichlet boundary data.
- Sphere mode solves on a triangulated unit sphere and compares against a known analytic eigenmode.
- The family API is surface-selected at comptime: one public module, two internal runtimes.
- Plane and sphere share a diffusion-specific time-stepped scaffold; mesh construction, exact solutions, and solve details stay local to each runtime.

## Verification

- Plane tests cover boundary pinning and second-order spatial convergence.
- Sphere tests cover monotone refinement improvement and an absolute analytic error bound.
- Short end-to-end checks live in [`examples/acceptance.zig`](../acceptance.zig).
