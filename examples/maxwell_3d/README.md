# maxwell_3d

3D cavity resonance example for the discrete exterior calculus Maxwell system on a tetrahedral mesh.

## Physics

The example evolves:

- `E ∈ Ω¹` on primal edges
- `B ∈ Ω²` on primal faces

with the semi-discrete update

- `∂B/∂t = -dE`
- `∂E/∂t = ★⁻¹ d ★ B`

inside a perfect electric conductor cavity. The magnetic constraint is structural: if `B` starts in `im(d)`, then `dB = ddA = 0` exactly, so the example asserts the milestone invariant through the update path rather than hoping truncation error stays small.

## Run

```sh
zig build -Doptimize=ReleaseFast example-maxwell3d -- --steps 1000 --dt 0.0025
zig build -Doptimize=ReleaseFast example-maxwell3d -- --steps 400 --dt 0.0025 --output output/maxwell3d --output-interval 40
```

The second command writes `.vtu` snapshots plus a `.pvd` collection file for ParaView.
