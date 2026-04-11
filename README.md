# flux

A type-safe PDE solver framework in Zig, built around discrete operators on
function spaces.

The core idea is simple: a solver is assembled by composing operators between
typed discrete fields on a mesh. flux provides the spaces, operators, and
composition rules; the user provides the graph.

See [project/vision.md](project/vision.md) for the full design philosophy.

| TE10 cavity resonance | Diffusion on a sphere |
|---|---|
| ![TE10 cavity animation](assets/cavity-512-grid-10000-steps.png) | ![Surface diffusion animation](assets/diffusion-sphere-r5-800-steps.png) |

## Status

Early development. The codebase is being shaped around exact structural
invariants and a small, generic public surface rather than backward-compatible
transitional APIs.

## Build

Requires Zig `0.15.2+`.

```sh
zig build
zig build test --summary all
zig build ci --summary all
zig build docs
```

## Run examples

The runnable physics surface lives in [examples/](examples/README.md).

```sh
zig build run -- maxwell --dim 2 --demo cavity
zig build run -- euler --dim 3 --steps 100
zig build run -- diffusion --surface sphere --refinement 3
```

Convenience run steps exist for each family:

```sh
zig build run-maxwell -- --dim 2 --demo dipole
zig build run-euler -- --dim 2 --demo dipole
zig build run-diffusion -- --surface plane --grid 32
```

Use `-Doptimize=ReleaseFast` for meaningful runs.

## Example families

| Family | Selector | Verification focus |
|---|---|---|
| [maxwell](examples/maxwell/README.md) | `--dim 2|3` | `dB = 0`, cavity convergence |
| [euler](examples/euler/README.md) | `--dim 2|3` | circulation and helicity preservation |
| [diffusion](examples/diffusion/README.md) | `--surface plane|sphere` | analytic convergence on flat and curved domains |

## Design commitments

- Geometric invariants are enforced structurally when possible.
- `comptime` carries mesh/form compatibility into the type system.
- All allocations are explicit.
- Tests are proof obligations, not smoke checks.

The exact invariants currently treated as first-class are:

- `dd = 0`
- `∇·B = 0`
- circulation conservation for incompressible Euler

## Development

The project is organized around epoch roadmaps, milestone acceptance criteria,
and issue-sized capabilities. The current architecture and process rules live
in:

- [project/vision.md](project/vision.md)
- [project/horizons.md](project/horizons.md)
- [project/components.md](project/components.md)

## Documentation

```sh
zig build docs
zig build serve-docs
```

## Contributing

The project is still in a heavy design/refactor phase and is not yet open to
outside contributions.
