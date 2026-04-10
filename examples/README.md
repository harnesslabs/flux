# Examples

Simulation examples demonstrating the flux operator stack across Maxwell,
Euler, and diffusion problems.

![TE10 cavity animation](../assets/cavity-512-grid-10000-steps.png)

---

## Quick start

For realistic throughput numbers, run the example in `ReleaseFast`:

```sh
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --demo dipole
zig build -Doptimize=ReleaseFast run-maxwell -- --dim 3 --steps 400 --dt 0.0025
zig build -Doptimize=ReleaseFast run-euler -- --dim 3 --steps 1000
zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --grid 32 --frames 8
```

The default `zig build` mode is `Debug`, which is useful for development but
materially slower on larger grids.

Visualize the output (requires Python 3.10+ and [uv](https://github.com/astral-sh/uv)):

```sh
uv run tools/visualize.py output --field B_flux --output animation.png
uv run tools/visualize.py output --field B_flux --output animation.gif
uv run tools/visualize.py output/euler_dipole --output animation.png
```

APNG is now the default recommendation because it preserves full color. Use
`.gif` only when you specifically need GIF compatibility; GIF is limited to a
256-color palette and will show more banding.

---

## Demos

| Demo | What it does | Example |
|------|-------------|---------|
| [**maxwell**](maxwell/README.md) | Discrete Maxwell examples in 2D and 3D with one command surface and shared field-evolution core | `zig build -Doptimize=ReleaseFast run-maxwell -- --dim 2 --steps 400 --frames 8` |
| [**euler**](euler/README.md) | Incompressible Euler examples in 2D and 3D with conserved circulation/helicity regression checks | `zig build -Doptimize=ReleaseFast run-euler -- --dim 2 --grid 32 --steps 1000` |
| [**diffusion**](diffusion/README.md) | Scalar diffusion on a plane or sphere with exact-solution convergence checks | `zig build -Doptimize=ReleaseFast run-diffusion -- --surface plane --grid 32 --frames 8` |

---

## CLI reference

```text
usage: flux <family> [options]

families:
  maxwell   use `--dim 2|3`
  euler     use `--dim 2|3`
  diffusion use `--surface plane|sphere`
```

### Parameter relationships

The timestep policy depends on the family:

- `maxwell --dim 2` derives `dt` from `--courant` and grid spacing unless
  `--dt` is passed explicitly.
- `euler --dim 2` derives `dt` from `--cfl` and grid spacing unless `--dt` is
  passed explicitly.
- `maxwell --dim 3` and `euler --dim 3` use explicit fixed `--dt`.
- `diffusion --surface plane` uses `dt_scale * h^2` unless `--dt` is passed.
- `diffusion --surface sphere` uses `final_time / steps`.

---

## Visualization tools

### tools/visualize.py — animated APNG/GIF

Parses VTK `.vtu` files and renders an animation with matplotlib. Triangle
meshes render as 2D pseudocolor plots; tetrahedral meshes render as 3D
cell-barycenter animations. No ParaView needed — `uv` handles dependencies
automatically.

```sh
uv run tools/visualize.py <input_dir> [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--field` | auto | Scalar field to plot; auto-selects `tracer`, `vorticity`, `B_flux`, `E_intensity`, or another available scalar |
| `--vectors` | auto | Optional 2D vector overlay; auto-selects `velocity` when present, use `none` to disable |
| `--output` | `<dir>/animation.png` | Output animation path; `.png`/`.apng` writes full-color APNG, `.gif` writes palette-limited GIF |
| `--fps` | 12 | Frames per second |

Available fields:

| Field | What it is |
|-------|-----------|
| `B_flux` | Magnetic flux per face (primal 2-form) |
| `E_intensity` | Electric field averaged to faces (primal 1-form projected) |
| `tracer` | Passive tracer transported by the 2D Euler showcase demo |
| `vorticity` | Face-centered 2D Euler vorticity |
| `stream_function` | Vertex-centered stream function |
| `velocity` | 2D Euler face velocity (vector overlay) |

### ParaView

For interactive exploration, open the `.pvd` collection file:

```sh
open output/cavity.pvd       # macOS
paraview output/dipole.pvd   # Linux
```

The `.pvd` indexes all snapshots with simulation timestamps for correct
time-axis playback. ParaView supports custom colormaps, streamlines, animation
export, and 3D viewing.

---

## Numerical background

The simulations use **Discrete Exterior Calculus (DEC)** on 2D simplicial
meshes. Fields are typed cochains — E is a primal 1-form (per edge), B is a
primal 2-form (per face) — with operator compatibility enforced at compile time.

The **leapfrog integrator** staggers E at integer steps and B at half-steps
(Yee-style), giving second-order time accuracy and symplectic structure (energy
conservation for source-free runs).

**Key invariant**: dB = 0 is enforced structurally — in 2D, B is a top-form, so
d₂B lives in the nonexistent 3-form space. The type system rejects
`exterior_derivative(B)` at compile time.

See the individual demo pages for physics details and parameter exploration
guides.
