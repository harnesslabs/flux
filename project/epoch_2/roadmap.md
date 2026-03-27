# Epoch 2: Library Architecture, N-Dimensional Generalization, and Physics Applications

Epoch 2 transforms flux from a working 2D Maxwell prototype into a general-purpose
FEEC/DEC library. The operator stack is hardened (convergent Hodge star, sparse
Laplacian, SIMD cochain arithmetic), generalized to arbitrary dimension (tested
concretely in 3D), and parameterized on Riemannian metrics. Physics-specific code
moves out of `src/` into `examples/` that consume flux as a Zig package — proving
the library architecture works by building five distinct simulations against it.

By the end of this epoch, flux can solve hyperbolic (Maxwell), parabolic (heat,
diffusion), and mixed elliptic-transport (incompressible Euler) problems on
imported meshes in 2D and 3D, with structural conservation guarantees verified
at every timestep.

---

## Milestones

### M1: Foundation Hardening + Library Architecture
**Goal:** Harden the operator stack for accuracy and performance, extract comptime
trait interfaces, establish flux as a consumable Zig package with physics code
living in `examples/`.
**Acceptance criterion:** `/audit` findings resolved. Whitney/Galerkin ★₁
convergence test passes at expected rate. `examples/maxwell_2d/` builds against
flux as a dependency and reproduces epoch 1 demos. Benchmarks establish baselines
with bounded CI checks. No physics-specific code remains in `src/`.
**Depends on:** none

Issues (target 8–10):
- Cochain carries mesh reference — cochains are typed functions on a complex, not bare float arrays
- Whitney/Galerkin mass matrix for ★₁ with convergence rate verification
- Sparse Laplacian assembly as symmetric sparse stiffness matrix
- Compressed incidence values for CsrMatrix — exploit {−1,0,+1} sparsity
- Comptime trait interfaces: Operator, TimeStepper, and Mesh concepts with compile-time contract checking
- Pluggable integrator extraction — leapfrog becomes a generic integrator parameterized on state type
- SIMD/vectorized batch operations on cochains (add, scale, negate, inner product)
- Benchmark suite with `zig build bench` and bounded CI regression checks
- Full `/audit` pass — safety, style, perf, tests — generate and resolve all findings
- Refactor `src/em/` → `examples/maxwell_2d/` consuming flux as a Zig package; library packaging in `build.zig`
- Close stale issues (#55, #39, #3)

### M2: N-Dimensional Operators + Mesh I/O + Linear Solver
**Goal:** Generalize all operators to arbitrary dimension (tested at n=2 and n=3),
add mesh import with auto-triangulation, parameterize the Hodge star on Riemannian
metrics, and implement a sparse linear solver.
**Acceptance criterion:** `dd = 0` and `★★⁻¹ = id` pass on random 3D tetrahedral
meshes. Poisson solve `Δu = f` converges at expected rate on an imported mesh.
All 2D tests pass unchanged (generic code, not duplicated). Riemannian ★ produces
correct results on a mesh with prescribed non-flat metric.
**Depends on:** M1 (trait interfaces, sparse Laplacian, library packaging)

Issues (target 8–10):
- Dimension-generic boundary operators (∂ₖ for arbitrary k via comptime, not separate implementations)
- Dimension-generic Hodge star (Whitney/Galerkin form generalized to n-simplices)
- Wedge product (∧) as a first-class operator with degree arithmetic checked at comptime
- Metric-parameterized Hodge star: `Metric(.flat)` (Euclidean) and `Metric(.riemannian)` (general tensor)
- Uniform tetrahedral grid constructor (3D analog of the 2D triangle grid)
- .obj mesh reader with automatic triangulation of quad faces
- .vtk/.vtu reader for unstructured grids
- 3D VTK export (tetrahedral cells + field data)
- Conjugate gradient solver for SPD systems with preconditioner interface (Jacobi as first implementation)
- Poisson solve: `Δu = f` with Dirichlet BCs on imported mesh, convergence test against analytic solution
- All property tests parameterized on dimension — dd=0, ★★⁻¹=id, ∂∂=0 tested at n=2 and n=3

### M3: Physics Applications
**Goal:** Five working simulations in `examples/`, exercising the full operator stack
across hyperbolic, parabolic, and mixed problems in 2D and 3D, each with structural
conservation guarantees.
**Acceptance criterion:** Discrete circulation conserved to machine precision over
1000 steps (Euler 2D/3D). ∇·B = 0 to machine precision at every timestep (Maxwell
3D). Heat equation converges at expected rate. Diffusion on curved surface matches
known analytic solution. All examples build as standalone binaries depending on
the flux library.
**Depends on:** M1 (library packaging, integrators), M2 (3D operators, solver, mesh I/O, Riemannian ★)

Issues (target 8–10):
- `examples/maxwell_3d/` — 3D cavity resonance on tetrahedral mesh, leapfrog integrator, ∇·B = 0 structural
- `examples/heat/` — heat equation via backward Euler + CG solve, convergence test against analytic solution
- `examples/euler_2d/` — vorticity-stream formulation, Poisson solve for stream function, circulation conservation on imported 2D mesh
- `examples/euler_3d/` — 3D vorticity via wedge product, helicity conservation as structural invariant
- `examples/diffusion_surface/` — heat equation on curved surface (sphere/torus) using Riemannian metric-dependent ★
- Boundary condition library as composable comptime types (PEC, no-slip, Dirichlet, periodic)
- Composable diagnostic observers (energy, circulation, helicity, divergence) — pluggable into any integrator
- Each example ships with: CLI interface, VTK output, convergence test, README explaining the physics and mathematics
