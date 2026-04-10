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
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo cavity --frames 50
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --steps 400 --dt 0.0025 --output output/maxwell3d --frames 10
```

## Notes

- The 2D mode supports `--demo dipole|cavity`, `--grid`, `--domain`,
  `--courant`, `--frequency`, and `--amplitude`.
- The 3D mode supports `--nx`, `--ny`, `--nz`, `--width`, `--height`,
  `--depth`, and `--output-interval`.
- The 3D tests verify exact preservation of `dB = 0` and convergence of the
  cavity reference mode under refinement.
