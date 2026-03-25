# TE₁₀ Cavity Resonance

A 2D electromagnetic standing wave in a square PEC (perfect electric conductor)
cavity — the simplest non-trivial test of the Maxwell solver against a known
analytical solution.

---

## Physics

A rectangular cavity with perfectly conducting walls supports discrete resonant
modes. The **TE₁₀ mode** (transverse electric, one half-wavelength in x, uniform
in y) on a `[0, L] × [0, L]` domain has:

| Field | Expression | DEC representation |
|-------|------------|--------------------|
| E_y(x,t) | sin(πx/L) · sin(ωt) | Primal 1-form (circulation on edges) |
| B_z(x,t) | cos(πx/L) · cos(ωt) | Primal 2-form (flux through faces) |
| E_x | 0 | — |

With wave speed c = 1:

- **Resonant frequency:** ω = π/L, so f = 1/(2L) Hz
- **Period:** T = 2L
- **Wavelength:** λ = 2L (one full wavelength fits the cavity)

The PEC boundary conditions (E tangential = 0 on walls) are automatically
satisfied: sin(πx/L) vanishes at x = 0 and x = L.

### What you see in the visualization

- **B_flux** (magnetic field): A cos(πx) pattern — positive (red) on the left
  half, negative (blue) on the right, with a nodal line at x = 0.5. This
  pattern oscillates in amplitude over time as energy transfers between E and B.

- **E_intensity** (electric field): A sin(πx) pattern — strongest at x = 0.5,
  zero at the walls. Quarter-period out of phase with B: when B is at maximum
  amplitude, E is zero, and vice versa.

The sawtooth texture visible in the plots comes from the triangular mesh — each
square cell is split by a SW→NE diagonal, creating a characteristic pattern in
the per-face field values.

### Energy conservation

Total electromagnetic energy U = ½(⟨E, ★₁E⟩ + ⟨B, ★₂B⟩) is conserved by the
symplectic leapfrog integrator. Over 2000 steps (~3 periods), energy drift is
< 0.02%. The energy oscillates between the E and B fields but the sum stays
constant — this is the defining property of a standing wave.

---

## Running the demos

### Prerequisites

- Zig 0.15.2+
- Python 3.10+ with [uv](https://github.com/astral-sh/uv) (for visualization)

### Cavity demo (source-free standing wave)

Initializes the exact TE₁₀ analytical mode and evolves it with the leapfrog
integrator. No source term — the standing wave persists indefinitely.

```sh
zig build run -- --demo cavity --steps 2000 --output output/cavity
```

### Dipole demo (default — driven simulation)

Places a sinusoidal point dipole at the cavity center, driving waves at the
TE₁₀ resonant frequency. Waves radiate outward and reflect off the PEC walls.

```sh
zig build run -- --steps 1000 --output output/dipole
```

### CLI reference

```
usage: flux [--demo <name>] [options]

demos:
  dipole    (default) point dipole radiating in a PEC cavity
  cavity    TE₁₀ standing wave in a PEC resonant cavity

options:
  --steps N     number of timesteps (default: 1000)
  --output DIR  output directory for VTK files (default: output)
  --help        show this message
```

Both demos use a 32×32 triangulated grid (2048 triangles) on a unit square
domain with dt = 0.003125 (Courant number 0.1). Output is a set of VTK `.vtu`
snapshot files plus a `.pvd` collection file for animation.

---

## Visualization

### Quick GIF (no ParaView needed)

```sh
uv run tools/visualize.py output/cavity --field B_flux --output cavity.gif
```

This parses the VTK XML files directly and renders each frame with matplotlib.
No ParaView installation required — `uv` handles all Python dependencies
automatically.

#### visualize.py reference

```
usage: uv run tools/visualize.py <input_dir> [options]

positional:
  input_dir            directory containing .vtu files

options:
  --field FIELD        CellData field to plot (default: B_flux)
  --output PATH        output GIF path (default: <input_dir>/animation.gif)
  --fps N              frames per second (default: 12)
```

Available fields in the VTK output:

| Field | Meaning |
|-------|---------|
| `B_flux` | Magnetic flux density (primal 2-form, per face) |
| `E_intensity` | Electric field intensity (primal 1-form projected to faces) |

### ParaView (full-featured)

If you have [ParaView](https://www.paraview.org/) installed, open the `.pvd`
collection file directly:

```sh
open output/cavity/cavity.pvd    # macOS
paraview output/cavity/cavity.pvd  # Linux
```

ParaView gives interactive 3D viewing, custom colormaps, streamlines, and
animation export. The `.pvd` file indexes all `.vtu` snapshots with their
simulation timestamps for correct time-axis playback.

---

## Numerical details

### Discretization

The simulation uses **Discrete Exterior Calculus (DEC)** on a 2D simplicial
mesh. Fields are represented as typed cochains — E is a primal 1-form
(one value per edge), B is a primal 2-form (one value per face). The type
system enforces correct operator composition at compile time.

The leapfrog time integrator staggers E at integer steps and B at half-steps
(Yee-style), giving second-order accuracy in time and exact conservation of the
symplectic structure.

### Eigenfrequency accuracy

The discrete eigenfrequency ω² converges to the analytical π²/L² with O(h²):

| Grid | h | ω²_numerical | ω²_analytical | Relative error |
|------|---|-------------|---------------|----------------|
| 8×8 | 0.125 | 9.838 | 9.870 | 0.32% |
| 16×16 | 0.0625 | 9.862 | 9.870 | 0.08% |

This is measured via the energy-based Rayleigh quotient ⟨d₁E, ★₂ d₁E⟩ / ⟨E, ★₁ E⟩,
which avoids the ★₁⁻¹ pseudo-inverse singularity on degenerate diagonal edges.

### Key invariant

The identity dB = 0 (no magnetic monopoles) is enforced **structurally**, not
numerically. In 2D, B is a top-form (2-form on a 2D mesh), so d₂B would live
in the 3-form space — which doesn't exist. The type system rejects
`exterior_derivative(B)` at compile time. This is the strongest possible
guarantee: the constraint cannot be violated because it cannot be expressed.
