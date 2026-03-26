# Dipole Radiation

A sinusoidal point source radiating electromagnetic waves inside a PEC cavity.
This is the default demo — the simplest way to see the Maxwell solver produce
non-trivial, dynamic field evolution.

```sh
zig build run -- --steps 1000 --output output/dipole
uv run tools/visualize.py output/dipole --field E_intensity --output dipole.gif
```

---

## Physics

A localized current source J oscillates sinusoidally on a single mesh edge near
the cavity center. This injects electromagnetic energy, which radiates outward
as waves, reflects off the PEC walls, and forms interference patterns.

The source is:

    J(t) = (A / edge_length) · sin(2πft)

concentrated on the edge nearest to the specified position (center of the
domain by default). All other edges have J = 0.

By default, the source frequency is set to the TE₁₀ resonant frequency
f = 1/(2L), which excites the cavity's fundamental mode. Changing `--frequency`
lets you drive at arbitrary frequencies — off-resonance driving produces more
complex interference patterns.

### What you see

- **E_intensity**: circular wavefronts expanding from the source point,
  reflecting off all four walls, and interfering. At the resonant frequency,
  the pattern eventually settles into something approximating the TE₁₀ standing
  wave. At off-resonance frequencies, it stays chaotic.

- **B_flux**: the magnetic field counterpart. Since B = −∫ dE dt, the B pattern
  lags E by a quarter cycle and shows complementary spatial structure.

---

## Parameters to explore

### Grid resolution

Finer grids resolve the wavefronts better and reduce numerical dispersion.
Compare the sharpness of the circular wavefronts.

```sh
# Coarse — fast, visible dispersion
zig build run -- --grid 16 --steps 500 --output output/dipole-16

# Fine — clean circular wavefronts
zig build run -- --grid 64 --steps 4000 --output output/dipole-64
```

### Source frequency

The default is the TE₁₀ resonance f = 1/(2L). Driving at other frequencies
changes the wavelength and interference pattern.

```sh
# Default: resonant frequency (f = 0.5 for L = 1)
zig build run -- --steps 2000 --output output/dipole-resonant

# Higher frequency — shorter wavelength, denser wavefronts
zig build run -- --frequency 2.0 --steps 2000 --output output/dipole-high-f

# Lower frequency — longer wavelength, broad pattern
zig build run -- --frequency 0.25 --steps 4000 --output output/dipole-low-f
```

### Amplitude

Controls the source strength. The physics is linear, so amplitude just scales
the field magnitudes. Useful for comparing energy growth rates.

```sh
zig build run -- --amplitude 0.1 --steps 1000 --output output/dipole-weak
zig build run -- --amplitude 10.0 --steps 1000 --output output/dipole-strong
```

### Domain size

A larger domain gives the wavefront more room to propagate before hitting walls.
The dt scales with the mesh spacing, so you may need more steps.

```sh
# Large domain — waves propagate further before reflecting
zig build run -- --domain 3.0 --grid 64 --steps 6000 --output output/dipole-large
```

### Courant number

Controls timestep size (dt = C · h). The default 0.1 is conservative. Larger
values run faster but introduce more numerical dispersion — the wavefronts
become less circular.

```sh
# See the effect of dispersion
zig build run -- --courant 0.4 --steps 500 --output output/dipole-dispersive
zig build run -- --courant 0.05 --steps 2000 --output output/dipole-accurate
```

---

## Energy growth

Unlike the cavity demo (which is source-free and conserves energy), the dipole
injects energy continuously. The energy printed at the end tells you how much
was injected. It should grow roughly as t² for on-resonance driving (resonant
buildup) and linearly for off-resonance.
