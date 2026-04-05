# Whitney Wedge Roadmap

This note records the intended direction after issue #84 lands.

## Current scope

Issue #84 adds a first-class primal-primal wedge product on cochains with the
API surface needed by the operator stack and the M3 vorticity work.

The implementation is intentionally scoped to the lowest-order FEEC/Whitney
setting:

- Interpret a cochain as coefficients in the Whitney basis.
- Take the smooth wedge in the interpolated space.
- Project back to a simplicial cochain by the de Rham map.
- Compute the induced local simplex formula directly instead of materializing
  the interpolated forms at runtime.

This keeps the user-facing operator simple while making the mathematical intent
explicit: the cochain product is induced from the Whitney/de Rham pipeline, not
invented ad hoc as a standalone array operation.

## Why this is not the end state

The true FEEC path is richer than "cochains in, cochains out".

The wedge of two Whitney forms is generally not itself a Whitney form. Projecting
back to cochains is the right first move for the current codebase, but it throws
away information that a higher-order FEEC representation could retain. That loss
shows up algebraically: graded commutativity and Leibniz survive, while strict
associativity does not in general.

For this reason, the current wedge should be understood as:

- the correct lowest-order operator for today's cochain stack;
- not the final nonlinear product story for flux.

## Recommended follow-on work

### 1. Introduce an explicit Whitney/FEEC form layer

Add types representing finite-element differential forms independently of the
cochain storage vectors. The key distinction is semantic:

- `Cochain` is discrete data attached to simplices.
- `WhitneyForm` or `FiniteElementForm` is a field representation with basis and
  evaluation semantics.

This would let the operator stack express "interpolate, operate, project"
without overloading `Cochain` as both storage and semantics.

### 2. Separate interpolation and projection operators

Make the two maps explicit:

- Whitney interpolation `W : C^k -> Λ_h^k`
- de Rham projection `R : Λ^k -> C^k`

Tracked by:
- #128 — explicit Whitney interpolation and de Rham projection operators

Today the wedge implementation effectively computes `R(W α ∧ W β)` directly.
That is fine for performance, but it hides the structure that future FEEC work
will need for higher-order spaces, metric variants, and more general nonlinear
operators.

### 3. Add higher-order FEEC spaces

Introduce polynomial differential-form families beyond the lowest-order Whitney
space. Once that exists, the project can represent nonlinear products without
immediate projection back to the lowest-order cochain complex.

This is the point where the wedge becomes a true FEEC operator instead of only
its lowest-order projected image.

Tracked by:
- #129 — finite-element form layer distinct from `Cochain`
- #130 — higher-order FEEC spaces for nonlinear operators

### 4. Generalize nonlinear operators together

The wedge should not evolve in isolation. The same interpolation/projection
infrastructure will likely be needed for:

- interior product,
- Lie derivative,
- nonlinear constitutive laws,
- vorticity/helicity diagnostics in M3 and beyond.

Treat these as one family of FEEC nonlinear operators rather than solving each
from scratch.

### 5. Add primal-dual and dual-dual wedge variants only after the FEEC layer exists

The current issue covers only the primal-primal product. Extending wedge to
mixed primal/dual inputs before the interpolation/projection story is explicit
would likely lock in the wrong abstraction boundary.

## Horizon alignment

This direction fits the existing project horizons:

- It keeps operators composable rather than baking nonlinear physics into one
  example.
- It stays compatible with generic scalar types and dimensionful quantities,
  because the interpolation/projection story is type-driven rather than tied to
  raw `f64` arrays.
- It avoids freezing the cochain layer as the only semantic layer, which is
  important if adaptive methods or richer FE spaces are added later.
