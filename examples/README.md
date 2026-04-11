# Examples

`examples/` is the runnable surface for the project. Each physics family has
one generic module and one CLI entry surface:

| Family | Selector | Notes |
|---|---|---|
| [maxwell](maxwell/README.md) | `--dim 2|3` | 2D dipole / cavity, 3D cavity reference mode |
| [euler](euler/README.md) | `--dim 2|3` | 2D circulation demo, 3D helicity regression mode |
| [diffusion](diffusion/README.md) | `--surface plane|sphere` | plane heat solve and spherical eigenmode |

## Quick start

Use `ReleaseFast` for real runs:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole --frames 8
zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --demo dipole --frames 8
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --grid 32 --frames 8
```

Use the umbrella binary when you want one command surface:

```sh
zig build run -- maxwell --dim 2 --demo cavity
zig build run -- euler --dim 3 --steps 100
zig build run -- diffusion --surface sphere --refinement 3
```

Ask a family for its own help:

```sh
zig build run -- maxwell --dim 2 --help
zig build run -- euler --dim 3 --help
zig build run -- diffusion --help
```

## Output

Examples write VTK snapshots and a `.pvd` series file into the chosen output
directory. Open the `.pvd` in ParaView or render an animation with:

```sh
uv run tools/visualize.py output --output animation.png
```

APNG is the default recommendation. Use `.gif` only when you specifically need
GIF compatibility.

## Shared CLI shape

- `--steps`, `--dt`, `--output`, and `--frames` are shared across families.
- `--grid` and `--domain` apply where the example is structured-grid based.
- `--refinement` and `--final-time` apply to the spherical diffusion path.

The timestep policy still depends on the family:

- Maxwell 2D uses `dt = courant * h` unless `--dt` is passed.
- Euler 2D uses `dt = cfl * h` unless `--dt` is passed.
- Diffusion plane uses `dt_scale * h^2` unless `--dt` is passed.
- Maxwell 3D, Euler 3D, and diffusion sphere use explicit fixed-time stepping.

## Verification

Each family keeps its own convergence/regression tests next to the example
module. [acceptance.zig](/Users/autoparallel/Code/flux/examples/acceptance.zig)
contains the short end-to-end milestone checks that run all invariants together.
