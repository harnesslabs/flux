# Epoch 1 Retrospective

Date: 2026-03-26
Duration: ~5 days (March 22 – March 26, 2026)

## What was planned

Three milestones building from mesh topology through operators to a working Maxwell simulation:

- **M1: Mesh + Visualization** (14 issues) — SoA mesh, boundary operators, geometry, VTK export
- **M2: Typed Forms + Discrete Operators** (14 issues) — comptime-typed cochains, d, ★, Laplacian
- **M3: Maxwell Simulation** (14 issues) — leapfrog FDTD, PEC BCs, dipole source, demos, CLI

Total: 42 issues, 3 milestones.

## What was delivered

Everything shipped. 42/42 issues closed, 20 PRs merged, 5,151 lines of Zig, 104 tests passing.

The framework goes from nothing to a working 2D Maxwell simulation with:
- Dimension-agnostic mesh with typed boundary operators
- Comptime-enforced k-form degree and primal/dual duality
- Exterior derivative, Hodge star, Laplace-de Rham operator
- VTK export with time-series support
- Leapfrog integrator with structural ∇·B = 0
- Two runnable demos (dipole radiation, cavity resonance) with CLI
- Python visualizer for GIF output

## Acceptance criteria status

| Milestone | Criterion | Status |
|-----------|-----------|--------|
| M1: Mesh + Visualization | ∂∂ = 0 exact; .vtu round-trips in ParaView | ✅ |
| M2: Typed Forms + Discrete Operators | Compile-time degree rejection; dd = 0 on 1000 random inputs; ★★⁻¹ = id | ✅ |
| M3: Maxwell Simulation | d₂B = 0 every timestep; 1000 steps stable; .vtu snapshots show field evolution | ✅ |

## What held

**Comptime type safety.** Degree and duality encoded as comptime parameters caught real bugs at compile time. Passing a dual cochain to a primal operator is a compile error. This is the right foundation — it should be extended, not relaxed.

**Property-based testing.** dd = 0 and ★★⁻¹ = id on 1000 random inputs caught issues that example-based tests would have missed. The fuzz tests are load-bearing correctness infrastructure.

**VTK pipeline.** Zero-dependency .vtu serializer with time-series support worked end-to-end. ParaView visualization was immediately useful for debugging the dipole radiation pattern.

**Standalone update functions.** Keeping `faraday_step` and `ampere_step` as free functions made it easy to compose the leapfrog integrator and will make pluggable integrators straightforward.

**The comptime-recursive `chain()`.** Type threading via `@TypeOf` on undefined inputs was elegant and caught real composition errors. No extra machinery needed.

## What didn't hold

**Decision logging — completely absent.** Zero entries during the epoch despite 6 non-trivial architectural choices. Every decision had to be reconstructed retroactively from PR descriptions. The `/decide` skill was never invoked. This is the biggest process failure of the epoch.

**Circumcentric dual assumption.** The initial choice of circumcentric dual for uniform grids propagated degeneracies (zero dual lengths on diagonal edges) through the Hodge star, Ampere step, and into the physics. Workarounds accumulated across 4+ PRs before the root cause was fixed in PR #69 by switching to barycentric dual. The degeneracy was documented in PR #47 as a "known limitation" rather than treated as a design problem to solve immediately. This violated the "zero technical debt" principle.

**Issue granularity was too fine.** Many issues were trivially batched: PR #47 closed 9 issues, PR #54 closed 5, PR #49 closed 3. Issues like "vertex coordinate storage" and "explicit allocator API" were not meaningful units of work — they were implementation details of the mesh struct. The roadmap called for 12–15 issues per milestone but several were below the "requires a non-obvious design choice" threshold.

**Milestones not closed.** All three milestones show 0 open issues but remain in "open" state on GitHub. Minor, but sloppy.

## Process observations

**Speed vs rigor tradeoff.** 42 issues in 5 days is fast. But the speed came at the cost of decision logging and issue quality. The epoch proved the architecture works end-to-end — that was the primary goal. But the process artifacts (decision log, issue sizing) need tightening for epoch 2.

**Workaround propagation.** The circumcentric dual issue is a case study in how a "known limitation" becomes technical debt. The fix was straightforward (swap circumcenter → barycenter) but wasn't done until the demo visibly broke. Earlier action would have avoided the pseudo-inverse workarounds in the Ampere step and the entry-by-entry degenerate-edge skipping in the Hodge star tests.

**PR batching was effective.** Bundling related issues into single PRs (e.g., all VTK export issues in PR #49, all Hodge star issues in PR #54) was the right call — splitting would have been churn. The issue templates should reflect that a PR can close multiple issues when they form a single logical unit.

**No `/decide` usage despite needing it.** The skill exists but was never triggered. Either the threshold for "non-obvious" was set too high in practice, or the friction of invoking `/decide` mid-implementation was too high.

## Recommendations for next epoch

1. **Enforce decision logging.** Lower the bar for what triggers `/decide`. Any new type, interface change, or workaround is a decision. Consider making `/decide` a checklist item on the PR template rather than relying on mid-implementation invocation.

2. **Coarser issue granularity.** Target 5–10 issues per milestone, not 12–15. Each issue should require at least one design choice and produce a reviewable PR. "Add field X to struct Y" is not an issue — it's a commit.

3. **Fix known limitations immediately.** "Known limitation" in a PR description should trigger either (a) a fix in the same PR, or (b) a blocking issue that prevents downstream work from building on the limitation. Do not document degeneracies and move on.

4. **Close milestones.** Automate or add to the PR-merge checklist: when the last issue in a milestone closes, close the milestone.

5. **Cochain should carry mesh reference.** Filed as #72. A cochain without its complex is just an array of floats — the mathematical semantics are lost. This should be addressed early in epoch 2.

6. **Whitney/Galerkin mass matrix for ★₁.** The barycentric diagonal Hodge star has ~8% non-converging error. Tracked as #70. This is the main accuracy bottleneck and should be a milestone-level goal.
