# Epoch 1: Foundations Through Maxwell Simulation

Epoch 1 builds flux from nothing to a working physics simulation. By the end, the framework
can represent a 2D triangulated mesh, define dimension-agnostic typed k-forms on it, apply
the discrete exterior derivative and Hodge star with exact invariant guarantees, export
results to VTK for ParaView inspection, and run a full Maxwell leapfrog simulation with
∇·B = 0 enforced structurally at every timestep. This is the proof that the architecture
works end-to-end — from topology through operators through physics through visualization.

The dimension design is 2D-first but n-agnostic: mesh and cochain types are parameterized
on dimension at the comptime level so 3D is a drop-in extension, not a rewrite.

---

## Milestones

### M1: Mesh + Visualization
**Goal:** Represent a 2D triangulated mesh with correct topology and export it — with
attached field data — to VTK format inspectable in ParaView.
**Acceptance criterion:** A uniform triangulated grid is constructed in memory; boundary
operators satisfy `∂∂ = 0` exactly on random meshes; a `.vtu` file is written with
0-form vertex data and 1-form edge data; the file parses in ParaView without errors and
values round-trip to machine precision.
**Depends on:** none

Issues (target 12–15):
- Mesh struct with SoA layout (`std.MultiArrayList` for vertices, edges, faces)
- Explicit allocator API (`init`, `deinit`, no hidden heap)
- Vertex coordinate storage (f64, n-dimensional via comptime)
- Oriented edge incidence: boundary operator ∂₁ in CSR format
- Oriented face incidence: boundary operator ∂₂ in CSR format
- Orientation consistency: verify ∂₁ ∘ ∂₂ = 0 exactly (property test on random meshes)
- Geometric data: edge lengths, face areas
- Dual cell volumes: circumcentric dual for 2D triangulations
- Uniform grid constructor: triangulated rectangle, configurable resolution
- VTK `.vtu` XML serializer (zero-dependency)
- 0-form → PointData export
- 1-form and 2-form → CellData export
- Time-series snapshot support (indexed `.vtu` files)
- Round-trip test: write then parse `.vtu`, verify values to machine precision

### M2: Typed Forms + Discrete Operators
**Goal:** Comptime-typed cochains with degree checked at compile time, exterior derivative
`d`, and Hodge star `★`, with mathematical invariants verified by property-based tests.
**Acceptance criterion:** A function accepting `Cochain(mesh, 1)` rejects `Cochain(mesh, 2)`
at compile time; `dd = 0` passes on 1000 random cochain inputs (fuzz-tested); `★ ∘ ★⁻¹ =
identity` holds to machine precision on random inputs.
**Depends on:** M1 (mesh geometry required for Hodge star volume ratios)

Issues (target 12–15):
- `Cochain(mesh, k)` type parameterized on mesh and degree via comptime
- Compile-time degree enforcement: incompatible compositions are compile errors
- Primal vs dual cochain distinction as comptime parameter
- Arithmetic on cochains: add, scale, negate
- Exterior derivative `d₀`: 0-form → 1-form (gradient)
- Exterior derivative `d₁`: 1-form → 2-form (curl)
- `dd = 0` fuzz property test: k=0 and k=1, 1000 random inputs
- Hodge star `★₀`: 0-form → n-form using dual vertex volumes / primal vertex volumes
- Hodge star `★₁`: 1-form → (n-1)-form using dual edge lengths / primal edge lengths
- Hodge star `★₂`: 2-form → 0-form (in 2D)
- Hodge star inverse `★⁻¹` (diagonal, exact)
- `★ ∘ ★⁻¹ = identity` property test on random inputs
- Laplace-Beltrami operator `Δ = d★d★ + ★d★d`
- Operator composition API: type-checked chaining of `d` and `★`
- Public API in `root.zig`: re-export `Cochain`, `d`, `star`, `laplacian`

### M3: Maxwell Simulation
**Goal:** End-to-end FDTD Maxwell simulation — radiating dipole on a 2D triangulated mesh
— with `d₂B = 0` verified at every timestep and ParaView-ready output.
**Acceptance criterion:** `d₂B = 0` holds to machine precision at every timestep; the
simulation runs 1000 steps without blowing up; `.vtu` snapshots are written every N steps
and open correctly in ParaView showing non-trivial field evolution.
**Depends on:** M1, M2

Issues (target 14–18):
- Field assignment: `E ∈ Ω¹` (primal 1-form on edges), `B ∈ Ω²` (primal 2-form on faces)
- Faraday update: `∂B/∂t = -dE` (half-step)
- Ampere-Maxwell update: `∂E/∂t = ★₁⁻¹ d ★₂ B - J`
- Leapfrog integrator: stagger E at integer steps, B at half-steps
- PEC boundary conditions: zero E on boundary edges
- Simulation state struct: `E`, `B`, `J`, timestep, mesh reference
- Current source `J`: point dipole, configurable frequency and position
- `d₂B = 0` invariant check: assert at every timestep (not just spot-check)
- Energy tracking: `‖E‖² + ‖B‖²` computed each step
- Simulation runner: configurable step count, output interval, output path
- VTK output integration: write E and B fields per snapshot
- TE/TM mode demo: 2D cavity resonance with known analytical modes
- Convergence test: halve grid spacing, verify field error decreases at expected rate
- CLI interface: `zig build run -- --steps N --output path/`
- End-to-end integration test: run 100 steps, assert ∇·B = 0, assert energy bounded
