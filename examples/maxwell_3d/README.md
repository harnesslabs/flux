# maxwell_3d

3D cavity resonance example for the discrete exterior calculus Maxwell system on a tetrahedral mesh.

## Physics

The example evolves:

- `E ∈ Ω¹` on primal edges
- `B ∈ Ω²` on primal faces

with the semi-discrete update

- `∂B/∂t = -dE`
- `∂E/∂t = ★⁻¹ d ★ B`

inside a perfect electric conductor cavity. The example initializes the analytical rectangular-cavity `TM110` mode through a discrete vector potential `A_z`, then sets `B = dA` on the mesh so the magnetic constraint is exact in the discrete complex as well as in the continuum.

The reference frequency is

```text
ω² = (π / width)² + (π / height)²
```

and the test suite checks that the discrete Rayleigh quotient converges toward this value under tetrahedral refinement.

## Run

```sh
zig build -Doptimize=ReleaseFast run-maxwell-3d -- --steps 1000 --dt 0.0025
zig build -Doptimize=ReleaseFast run-maxwell-3d -- --steps 400 --dt 0.0025 --output output/maxwell3d --output-interval 40
uv run tools/visualize.py output/maxwell3d --field B_flux --output output/maxwell3d/animation.png
```

The second command writes `.vtu` snapshots plus a `.pvd` collection file for ParaView. The third command uses the shared visualizer, which now renders tetrahedral output as a 3D barycenter animation.

## CLI reference

```text
usage:
  zig build -Doptimize=ReleaseFast run-maxwell-3d -- [options]

mesh:
  --nx N              tetrahedral cells in x (default: 2)
  --ny N              tetrahedral cells in y (default: 2)
  --nz N              tetrahedral cells in z (default: 2)
  --width L           cavity width  (default: 1.0)
  --height L          cavity height (default: 1.0)
  --depth L           cavity depth  (default: 1.0)

time stepping:
  --steps N           leapfrog steps (default: 1000)
  --dt DT             fixed timestep (default: 0.01)

output:
  --output DIR        write VTK snapshots into DIR
  --output-interval N write every N steps when output is enabled
```

The terminal output matches the 2D example style: run header, live progress bar, and a final results block with `||dB||₂`, elapsed time, and visualization hint.
