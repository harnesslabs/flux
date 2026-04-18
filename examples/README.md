# Examples

`examples/` is the runnable surface for the project. Each physics family has
one generic module and one CLI entry surface:

| Family | Selector | Notes |
|---|---|---|
| `new_maxwell` | `maxwell <scenario>` | `dipole` or `cavity`, with `cavity --dim 2|3` |
| `new_euler` | `euler <scenario>` | `gaussian`, `dipole`, or `reference` |
| `new_diffusion` | `diffusion <scenario>` | `plane` or `sphere` |

## Quick start

Use `ReleaseFast` for real runs:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --width 1.0 --height 1.0 --frames 8
zig build -Doptimize=ReleaseFast run-euler -- --width 1.0 --height 1.0 --frames 8
zig build -Doptimize=ReleaseFast run-diffusion -- --width 1.0 --height 1.0 --frames 8
```

Use the umbrella binary when you want one command surface:

```sh
zig build run -- maxwell cavity --dim 2
zig build run -- euler reference --steps 100
zig build run -- diffusion sphere --refinement 3
```

Ask a family for its own help:

```sh
zig build run -- maxwell cavity --help
zig build run -- euler gaussian --help
zig build run -- diffusion plane --help
```

## Output

Examples write VTK snapshots and a `.pvd` series file into the chosen output
directory. Open the `.pvd` in ParaView or render an animation with:

```sh
uv run tools/visualize.py output --output animation.png
```

APNG is the default recommendation. Use `.gif` only when you specifically need
GIF compatibility.

## CLI shape

- `maxwell dipole` is the default Maxwell convenience scenario behind `run-maxwell`.
- `euler gaussian` is the default Euler convenience scenario behind `run-euler`.
- `diffusion plane` is the default diffusion convenience scenario behind `run-diffusion`.
- `--steps`, `--dt`, `--output`, `--frames`, and `--output-interval` are the shared execution flags.
- Per-family geometry/resolution flags are scenario-specific by design rather than forced through one shared schema.

## Verification

Each family keeps its own convergence/regression tests next to the example
module. [acceptance.zig](/Users/autoparallel/Code/flux/examples/acceptance.zig)
contains the short end-to-end milestone checks that run all invariants together.
