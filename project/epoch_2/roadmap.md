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
**GitHub milestone:** https://github.com/harnesslabs/flux/milestone/4

Issues (10):
- [x] ~~#55 — Investigate barycentric dual~~ (closed: investigation complete, fix tracked by #70)
- [ ] #72 — Cochain carries mesh reference — cochains are typed functions on a complex, not bare float arrays
- [ ] #70 — Whitney/Galerkin mass matrix for ★₁ with convergence rate verification
- [ ] #58 — Sparse Laplacian assembly as symmetric sparse stiffness matrix
- [ ] #48 — Compressed incidence values for CsrMatrix — exploit {−1,0,+1} sparsity
- [ ] #62 — Comptime TimeStepper interface — enforce integrator contract at compile time
- [ ] #75 — Comptime Operator and Mesh trait interfaces with compile-time contract checking
- [ ] #76 — Extract pluggable integrator — leapfrog as generic time stepper parameterized on state type
- [ ] #77 — SIMD/vectorized batch operations on cochains (add, scale, negate, inner product)
- [ ] #78 — Benchmark suite with `zig build bench` and bounded CI regression checks
- [ ] #79 — Refactor `src/em/` → `examples/maxwell_2d/` consuming flux as a Zig package

Pre-milestone housekeeping:
- [ ] #3 — Add/fix branch protections (priority/high)
- #39 — d₂B = 0 structural assertion → deferred to M2 (meaningful in dim > 2)

### M2: N-Dimensional Operators + Mesh I/O + Linear Solver
**Goal:** Generalize all operators to arbitrary dimension (tested at n=2 and n=3),
add mesh import with auto-triangulation, parameterize the Hodge star on Riemannian
metrics, and implement a sparse linear solver.
**Acceptance criterion:** `dd = 0` and `★★⁻¹ = id` pass on random 3D tetrahedral
meshes. Poisson solve `Δu = f` converges at expected rate on an imported mesh.
All 2D tests pass unchanged (generic code, not duplicated). Riemannian ★ produces
correct results on a mesh with prescribed non-flat metric.
**Depends on:** M1 (trait interfaces, sparse Laplacian, library packaging)
**GitHub milestone:** https://github.com/harnesslabs/flux/milestone/5

Issues (10):
- [ ] #80 — Generalize Mesh topological dimension from constant 2 to comptime parameter
- [ ] #81 — Uniform tetrahedral grid constructor for 3D meshes
- [ ] #82 — Dimension-generic Hodge star with n-simplex geometric formulas
- [ ] #83 — Dimension-generic exterior derivative with 3D dd=0 property tests
- [ ] #84 — Wedge product (∧) as first-class operator with comptime degree arithmetic
- [ ] #85 — Metric-parameterized Hodge star — flat (Euclidean) and Riemannian modes
- [ ] #86 — .obj mesh reader with separable triangulation pass
- [ ] #87 — 3D VTK export for tetrahedral meshes with field data
- [ ] #88 — Conjugate gradient solver for SPD systems with preconditioner interface
- [ ] #89 — Poisson solve Δu = f with convergence rate verification

Design note: The .obj reader (#86) produces a `RawMesh` with polygon faces of
any valence, then an explicit `triangulate()` pass converts to simplicial form.
Native quad/polygon cell support is tracked as a horizon item (see `horizons.md`).

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

**GitHub milestone:** https://github.com/harnesslabs/flux/milestone/6

Issues (9):
- [ ] #90 — Composable comptime boundary condition types (PEC, Dirichlet, no-slip, periodic)
- [ ] #91 — Composable diagnostic observer framework pluggable into any integrator
- [ ] #92 — `examples/maxwell_3d/` — 3D cavity resonance with ∇·B = 0 structural guarantee
- [ ] #93 — `examples/heat/` — heat equation via backward Euler + CG solve with convergence verification
- [ ] #94 — `examples/euler_2d/` — vorticity-stream formulation with circulation conservation
- [ ] #95 — `examples/euler_3d/` — 3D vorticity via wedge product with helicity conservation
- [ ] #96 — `examples/diffusion_surface/` — heat equation on curved surface with Riemannian Hodge star
- [ ] #97 — Example packaging infrastructure — shared CLI, VTK output, and build integration
- [ ] #98 — End-to-end acceptance test verifying all M3 conservation invariants in CI
