# Epoch 1 Decision Log

<!-- Append entries with /decide. Format: ## YYYY-MM-DD: <title> -->

## 2026-03-22: SW-NE diagonal split for uniform grid [retroactive]

**Decision:** Uniform grid triangulation uses a consistent SW-NE diagonal split per cell, producing right triangles.

**Alternatives considered:**
- Alternating diagonal splits (would avoid systematic degeneracy)
- Delaunay-quality triangulation (more complex, better element quality)

**Rationale:** Simplest scheme. But right triangles sharing a hypotenuse have coincident circumcenters, which created the degenerate dual edge problem that propagated through ★₁, the Ampere step, and the dipole demo. In hindsight, the mesh generator should not have been the source of operator degeneracies.

**Source:** PR #47, commit a1cec65

## 2026-03-22: Circumcentric dual → barycentric dual [retroactive]

**Decision:** Switched `uniform_grid` from circumcentric to barycentric dual after the circumcentric dual produced `dual_length = 0` on all diagonal edges, making ★₁ singular on ~1/3 of edges.

**Alternatives considered:**
- Keep circumcentric dual with pseudo-inverse workarounds (chosen initially, propagated workarounds into Hodge star, Ampere step, and tests)
- Non-degenerate mesh generator (would fix root cause but doesn't help uniform grids)
- Whitney/Galerkin mass matrix (correct long-term fix, tracked as #70)

**Rationale:** Should have chosen barycentric from the start. The circumcentric degeneracy was a known issue from PR #47 but was worked around rather than fixed. The workarounds polluted multiple modules (Hodge star round-trip tests, Ampere step pseudo-inverse logic). Barycentric dual trades primal-dual orthogonality for non-degeneracy — introduces ~8% systematic error that doesn't converge, but this is acceptable until the Whitney/Galerkin mass matrix replaces the diagonal ★₁.

**Source:** PR #69, commit 0a44562

## 2026-03-22: Cochain carries mesh pointer [retroactive, corrected 2026-03-28]

**Decision:** `Cochain` stores `values: []f64` and `mesh: *const MeshType`. Operators extract the mesh from the cochain directly — no separate mesh argument.

**Alternatives considered:**
- Pure data cochain with no mesh reference (considered in PR #51 review, rejected)

**Rationale:** A cochain is mathematically a function on the cells of a CW complex — it is inherently tied to its topology. Carrying the mesh reference is faithful to the mathematics and eliminates boilerplate at operator call sites. The mesh pointer was retained in PR #51 (commit e51d0e1: "keep mesh pointer"). The original decision log entry incorrectly stated the pointer was removed.

**Source:** PR #50, PR #51 (commit e51d0e1, fc10c7b)

## 2026-03-24: Standalone update functions, not methods on State [retroactive]

**Decision:** `faraday_step` and `ampere_step` are free functions, not methods on `MaxwellState`.

**Alternatives considered:**
- Method-style `state.faraday_step(dt)`

**Rationale:** Separates state representation from integration logic, aligning with the pluggable time-integrator horizon. Reasonable for now, but as more physics modules appear, a generic `State` type with methods may be preferable to ad-hoc free functions per physics domain.

**Source:** PR #60, commit 81e2416

## 2026-03-24: Direct diagonal Ampere step with pseudo-inverse [retroactive]

**Decision:** Ampere step used manual diagonal ★₁⁻¹ with pseudo-inverse (0/0 → 0) instead of the generic `compose.chain` with `hodge_star_inverse`.

**Alternatives considered:**
- Use generic chain composition (panics on zero entries)

**Rationale:** Workaround for the circumcentric dual degeneracy. Superseded by PR #69's barycentric dual switch, which eliminated the zero entries entirely. This decision was a symptom of the earlier wrong choice on circumcentric dual.

**Source:** PR #60, commit 81e2416; removed in PR #69

## 2026-03-24: Comptime-recursive type threading for chain() [retroactive]

**Decision:** `chain` uses `@TypeOf(op(undefined_allocator, undefined_input))` to resolve each operator's return type at comptime without evaluating function bodies.

**Alternatives considered:**
- Trait/interface-based type registration
- Runtime type checks

**Rationale:** Leverages Zig's type system naturally — degree/duality mismatches become compile errors through operator signatures alone. Good choice for now. Trait/interface abstractions may be added later as the operator algebra grows, but the comptime approach remains the right foundation.

**Source:** PR #59, commit d703bda

## 2026-03-28: Whitney assembly maps local edges to global edges by vertex lookup

**Decision:** The boundary_2 column → local edge mapping in Whitney mass matrix assembly is built by explicit vertex-pair matching, not by assuming column order matches local edge order.

**Alternatives considered:**
- Assume boundary_2 columns follow the face's local edge order (fast, but wrong — boundary_2 stores edges sorted by global index)

**Rationale:** The boundary_2 row stores edges sorted by global edge index, which coincidentally matches the lower-right triangle's local edge order but not the upper-left triangle's (where the diagonal edge is local edge 0 but appears last in boundary_2). Matching by vertex pair is O(9) per face and costs nothing compared to assembly. Correctness over cleverness.

**Source:** PR #105, issue #70
