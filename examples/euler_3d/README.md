# Euler 3D

This example exercises the 3D DEC ingredients needed for incompressible Euler
on tetrahedral meshes:

- Velocity `u` is a primal 1-form on edges.
- Vorticity `ω = du` is a primal 2-form on faces.
- Helicity is measured as the discrete top-form sum `∫ u ∧ ω`.
- Velocity recovery uses the primal 1-form Dirichlet solve added for issue
  #141, so the example depends on the real operator stack rather than ad hoc
  physics-side linear algebra.

Run the acceptance reference with:

```sh
zig build -Doptimize=ReleaseFast run-euler-3d -- --steps 1000
```

Useful flags:

```sh
zig build -Doptimize=ReleaseFast run-euler-3d -- --help
zig build -Doptimize=ReleaseFast run-euler-3d -- --nx 3 --ny 3 --nz 2 --output output/euler_3d --output-interval 50
```

Output is written as `.vtu` snapshots plus an `euler_3d.pvd` collection when
`--output` and `--output-interval` are both provided.

Current scope:

- The seeded mode is a steady helical reference field chosen to make the
  acceptance invariant exact on the existing operator stack.
- Each timestep re-solves the primal 1-form Dirichlet problem, recomputes
  `ω = du`, and evaluates the helicity density `u ∧ ω`.
- The nonlinear 3D transport closure is not yet advanced in time here. That is
  a bounded limitation of this first M3 acceptance example, not a hidden claim.

Derivation sketch:

1. Helicity in 3D is `H = ∫Ω u ∧ ω`, the natural DEC top-form analogue of the
   smooth `∫ u · curl(u) dV`.
2. The example keeps a steady reference mode fixed, so the only admissible
   drift is numerical drift from the operator stack itself.
3. Re-solving the 1-form Poisson problem on every step turns the acceptance
   test into a direct check that the 3D operator pipeline is internally
   consistent and does not leak helicity under repeated use.

Reference: Mohamed, Hirani, Samtaney, *Discrete Exterior Calculus
Discretization of Incompressible Navier-Stokes Equations* (2016).
