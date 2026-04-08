# Euler 2D

This example advances incompressible 2D Euler on a triangular mesh in a
vorticity-stream formulation:

- Vorticity `ω` is stored per triangle.
- The stream function `ψ` is recovered from `Δψ = ω` with the library Poisson solve.
- Face velocity is reconstructed as the rotated gradient `u = ★dψ`.
- Vorticity is transported by a conservative face-flux update, so total
  circulation `∫Ω ω dA` is preserved by pairwise internal flux cancellation.

Run the invariant reference with:

```sh
zig build -Doptimize=ReleaseFast run-euler-2d -- --demo gaussian --grid 32 --steps 1000
```

For an actually dynamic showcase animation, use the vortex dipole demo:

```sh
zig build -Doptimize=ReleaseFast run-euler-2d -- --demo dipole --grid 48 --steps 1200
uv run tools/visualize.py output/euler_dipole
```

Useful flags:

```sh
zig build -Doptimize=ReleaseFast run-euler-2d -- --help
zig build -Doptimize=ReleaseFast run-euler-2d -- --frames 0
```

Output is written as `.vtu` snapshots plus an `euler_2d.pvd` collection in
`output/euler_2d/` by default for the Gaussian demo and `output/euler_dipole/`
for the dipole demo. Open the `.pvd` file in ParaView or use the existing
`tools/visualize.py` helper on the snapshot directory.

Demo modes:

- `gaussian` keeps the current near-steady single-vortex setup for the
  circulation-invariant reference case.
- `dipole` initializes a translating vortex pair and a passive tracer stripe.
  The tracer is advected by the same face-flux transport step, which makes the
  motion visually obvious in the rendered animation.

Derivation sketch:

1. In 2D incompressible flow, velocity is represented through a scalar stream
   function `ψ`, with `u = ∇⊥ψ`.
2. Vorticity is `ω = ∇×u = Δψ`, so recovering `ψ` each timestep reduces to a
   Poisson solve with homogeneous Dirichlet boundary data.
3. The transport equation `∂tω + u·∇ω = 0` is discretized here as a
   conservative upwind flux between neighboring triangles. That discretization
   is only first order in space, but it enforces the issue's key invariant:
   the sum of face circulation masses changes only by boundary flux, which is
   zero in this closed domain.

Reference: Elcott et al., *Stable, Circulation-Preserving, Simplicial Fluids*
(2007).
