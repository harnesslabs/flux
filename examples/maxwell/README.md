# Maxwell

Discrete exterior calculus Maxwell examples with one CLI surface:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --help
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --help
```

## Dimensions

- `--dim 2` runs the structured simplicial cavity and dipole examples.
- `--dim 3` runs the tetrahedral rectangular-cavity `TM110` reference mode.

Both dimensions evolve:

- `E ∈ Ω¹` on primal edges
- `B ∈ Ω²` on primal faces

with the semi-discrete update

- `∂B/∂t = -dE`
- `∂E/∂t = ★⁻¹ d ★ B - J`

The shared leapfrog/Faraday/boundary plumbing lives in the common Maxwell
module so the command surface and top-level build graph no longer split by
dimension.

## Run

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole --grid 32 --steps 800
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo cavity --frames 50
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --steps 400 --dt 0.0025 --output output/maxwell3d --frames 10
```

## 2D demos

### `dipole`

- Drives one localized current source near the cavity center.
- Good for visual sanity checks: outgoing waves, wall reflections, and
  resonance buildup.
- Main tuning knobs are `--frequency`, `--amplitude`, `--grid`, and
  `--courant`.

Typical runs:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole --grid 16 --steps 500
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole --frequency 2.0 --steps 1000
```

### `cavity`

- Seeds the TE10 standing-wave mode directly and runs source-free.
- This is the clean regression path for energy drift and eigenvalue
  convergence.
- Main tuning knobs are `--grid`, `--domain`, and `--courant`.

Typical runs:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo cavity --grid 8 --steps 500
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo cavity --grid 64 --steps 4000
```

## 3D mode

- Seeds the rectangular-cavity `TM110` reference mode on a tetrahedral mesh.
- This is a verification path, not a menu of separate demos.
- Main tuning knobs are `--nx`, `--ny`, `--nz`, `--dt`, and `--output-interval`.

Typical run:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --nx 4 --ny 4 --nz 4 --steps 400 --dt 0.0025
```

## Notes

- The 2D mode supports `--demo dipole|cavity`, `--grid`, `--domain`,
  `--courant`, `--frequency`, and `--amplitude`.
- The 3D mode supports `--nx`, `--ny`, `--nz`, `--width`, `--height`,
  `--depth`, and `--output-interval`.
- Convergence tests live beside the family module and cover TE10 energy/eigenvalue
  behavior in 2D and TM110 refinement in 3D.
