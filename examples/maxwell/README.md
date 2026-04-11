# Maxwell

One family module, one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --help
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --help
```

## Modes

| Selector | Purpose | Main knobs |
|---|---|---|
| `--dim 2 --demo dipole` | driven cavity run for visual checks | `--grid`, `--courant`, `--frequency`, `--amplitude` |
| `--dim 2 --demo cavity` | TE10 standing-wave regression path | `--grid`, `--domain`, `--courant` |
| `--dim 3` | TM110 rectangular-cavity reference mode | `--nx --ny --nz`, `--dt`, `--output-interval` |

## Run

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole --grid 32 --steps 800
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo cavity --grid 32 --steps 2000
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --nx 4 --ny 4 --nz 4 --steps 400 --dt 0.0025
```

## Model

- Fields are primal cochains: `E ∈ Ω¹`, `B ∈ Ω²`.
- The semi-discrete update is `∂B/∂t = -dE`, `∂E/∂t = ★⁻¹ d ★ B - J`.
- 2D and 3D share the same family API even though their internal runtimes differ.

## Verification

- 2D convergence tests cover TE10 eigenvalue accuracy and bounded energy drift.
- 3D convergence tests cover TM110 refinement behavior.
- Short end-to-end invariant checks live in [`examples/acceptance.zig`](../acceptance.zig).
