# Euler

Incompressible Euler examples with one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --help
zig build -Doptimize=ReleaseFast run-euler -- --dim 3 --help
```

## Dimensions

- `--dim 2` runs the vorticity-stream example on a triangular mesh.
- `--dim 3` runs the tetrahedral helical reference mode used for helicity
  regression checks.

## 2D

- Vorticity `ω` is stored per triangle.
- The stream function `ψ` is recovered from `Δψ = ω`.
- Face velocity is reconstructed as the rotated gradient `u = ★dψ`.
- Conservative face-flux transport preserves total circulation.

Example runs:

```sh
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --demo gaussian --grid 32 --steps 1000
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --demo dipole --grid 48 --steps 1200
```

## 3D

- Velocity `u` is a primal 1-form on edges.
- Vorticity `ω = du` is a primal 2-form on faces.
- Helicity is measured as `∫ u ∧ ω`.
- Velocity recovery uses the library 1-form Dirichlet solve.

Example run:

```sh
zig build -Doptimize=ReleaseFast run-euler -- --dim 3 --steps 1000
```

The seeded 3D mode is intentionally steady. It is a regression harness for the
current operator stack, not a claim that the full nonlinear 3D closure is
already advanced here.
