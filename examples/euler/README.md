# Euler

One family module, one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --help
zig build -Doptimize=ReleaseFast run-euler -- --dim 3 --help
```

## Modes

| Selector | Purpose | Main knobs |
|---|---|---|
| `--dim 2 --demo gaussian` | smooth vorticity-stream transport demo | `--grid`, `--cfl`, `--steps` |
| `--dim 2 --demo dipole` | stronger circulation-preserving transport case | `--grid`, `--cfl`, `--steps` |
| `--dim 3` | steady helical reference mode for helicity regression | `--nx --ny --nz`, `--dt`, `--output-interval` |

## Run

```sh
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --demo gaussian --grid 32 --steps 1000
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --demo dipole --grid 48 --steps 1200
zig build -Doptimize=ReleaseFast run-euler -- --dim 3 --steps 1000
```

## Model

- 2D stores vorticity per face, recovers a stream function, and reconstructs face velocity as the rotated gradient.
- 3D stores velocity as a primal 1-form, vorticity as `du`, and measures helicity `∫ u ∧ ω`.
- The 3D seeded mode is intentionally steady. It is a regression harness for the current operator stack, not a full nonlinear 3D evolution claim.

## Verification

- 2D example runs preserve total circulation.
- 3D example runs preserve helicity for the seeded reference mode.
- Short end-to-end invariant checks live in [`examples/acceptance.zig`](../acceptance.zig).
