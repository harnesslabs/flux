# TE₁₀ Cavity Resonance

A source-free electromagnetic standing wave in a square PEC cavity — the
simplest non-trivial validation of the Maxwell solver against a known
analytical solution.

```sh
zig build run -- --demo cavity --steps 2000 --output output/cavity
uv run tools/visualize.py output/cavity --field B_flux --output cavity.gif
```

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

- **Resonant frequency:** ω = π/L, giving f = 1/(2L) Hz
- **Period:** T = 2L
- **Wavelength:** λ = 2L (one full wavelength fits the cavity)

The PEC boundary conditions (tangential E = 0 on walls) are automatically
satisfied: sin(πx/L) vanishes at x = 0 and x = L.

### Initial conditions

The demo initializes the exact analytical mode and runs **source-free** — no
current drives the system. The leapfrog integrator uses staggered fields:

- **E** at t = 0: sin(πx/L) · sin(0) = **0** (zero-initialized)
- **B** at t = −dt/2: cos(πx/L) · cos(−ωdt/2) ≈ cos(πx/L) (projected onto face integrals)

The standing wave persists indefinitely, trading energy back and forth between
E and B.

### What you see

- **B_flux**: cos(πx/L) pattern — red on the left, blue on the right, nodal
  line at x = L/2. Oscillates in amplitude as energy flows between E and B.
- **E_intensity**: sin(πx/L) pattern — strongest at x = L/2, zero at walls.
  Quarter-period out of phase with B.

The sawtooth texture in the plots comes from the triangulated mesh — each square
cell is split by a SW→NE diagonal.

---

## Parameters to explore

### Grid resolution: coarse vs fine

Compare the mesh-induced artifact pattern at different resolutions. The standing
wave is the same; the triangulation effects shrink with refinement.

```sh
# Coarse (visible triangulation artifacts)
zig build run -- --demo cavity --grid 8 --steps 500 --output output/cavity-8x8

# Medium (default)
zig build run -- --demo cavity --grid 32 --steps 2000 --output output/cavity-32x32

# Fine (smooth, but slower)
zig build run -- --demo cavity --grid 64 --steps 4000 --output output/cavity-64x64
```

### Domain size

Changing the domain length L changes the resonant frequency (ω = π/L). A larger
cavity resonates at a lower frequency with a longer period.

```sh
# Small cavity — high frequency, fast oscillation
zig build run -- --demo cavity --domain 0.5 --steps 2000 --output output/cavity-small

# Large cavity — low frequency, slow oscillation
zig build run -- --demo cavity --domain 2.0 --steps 2000 --output output/cavity-large
```

### Courant number

Controls the timestep size relative to the mesh spacing (dt = C · h). Smaller
values are more accurate but need more steps to cover the same physical time.

```sh
# Conservative (very accurate, slow)
zig build run -- --demo cavity --courant 0.05 --steps 4000 --output output/cavity-conservative

# Aggressive (faster, more dispersion)
zig build run -- --demo cavity --courant 0.4 --steps 500 --output output/cavity-fast
```

---

## Energy conservation

Total electromagnetic energy U = ½(⟨E, ★₁E⟩ + ⟨B, ★₂B⟩) is conserved by the
symplectic leapfrog integrator. The energy ratio printed at the end should stay
close to 1.0 — any drift indicates the Courant number is too large for the mesh.

| Run | Energy ratio after 3 periods |
|-----|------------------------------|
| 32×32, C = 0.1 | 0.999818 |
| 64×64, C = 0.1 | 0.999349 (500 steps, partial period) |

---

## Eigenfrequency accuracy

The discrete eigenfrequency ω² converges to the analytical π²/L² with O(h²):

| Grid | h | ω²_numerical | ω²_analytical | Relative error |
|------|---|-------------|---------------|----------------|
| 8×8 | 0.125 | 9.838 | 9.870 | 0.32% |
| 16×16 | 0.0625 | 9.862 | 9.870 | 0.08% |

Measured via the energy-based Rayleigh quotient ⟨d₁E, ★₂ d₁E⟩ / ⟨E, ★₁ E⟩,
which avoids the ★₁⁻¹ pseudo-inverse singularity on degenerate diagonal edges.
This test runs in CI (`zig build test`).
