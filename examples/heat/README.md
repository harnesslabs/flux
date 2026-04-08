# Heat Example

Implicit heat-equation example on the unit square:

\[
\partial_t u = \Delta u,\qquad
u(x,y,0) = \sin(\pi x)\sin(\pi y),\qquad
u|_{\partial \Omega} = 0.
\]

The implementation uses backward Euler in time and solves the SPD system

\[
(M + \Delta t\,S)u^{n+1} = M u^n
\]

with conjugate gradient on the interior vertices, where `M` is the diagonal
DEC `★₀` mass and `S` is the assembled 0-form stiffness matrix.

## Run

```sh
zig build -Doptimize=ReleaseFast run-heat -- --grid 32 --frames 8
zig build -Doptimize=ReleaseFast run-heat -- --help
```

## Output

Snapshots are written to `output/heat/` by default:

- `heat_XXXX.vtu` stores point fields `temperature`, `temperature_exact`, and `temperature_error`
- `heat.pvd` indexes the time series for ParaView

## Verification

The example test suite verifies:

- zero initial data stays zero under backward Euler
- homogeneous Dirichlet boundary values remain pinned to zero
- the spatial error converges at second order on successively refined grids

The analytical reference is

\[
u(x,y,t)=e^{-2\pi^2 t}\sin(\pi x)\sin(\pi y).
\]
