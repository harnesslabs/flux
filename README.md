# flux

A composable, type-safe PDE solver framework in Zig.

The central abstraction is the **operator on a function space**: a map between two spaces of discrete fields defined on a mesh. A PDE solver is a directed acyclic graph of such operators. flux provides the nodes (function spaces), the edges (operators), and the composition rules. The user provides the graph.

See [`project/vision.md`](project/vision.md) for the full design philosophy.

## Status

Early development. See [Projects](https://github.com/harnesslabs/flux/projects) for the current epoch roadmap.

---

## Building

Requires Zig 0.15.2+.

```sh
zig build           # compile
zig build test      # run all tests
zig build run       # run the CLI
zig build docs      # generate API docs to zig-out/docs/
```

---

## Design

### Abstraction hierarchy

Function spaces and operators compose into solver graphs. Each level is independently swappable — changing a discretization or time integrator touches one node, nothing else.

### Core invariants

These hold exactly (to machine precision) and are verified by property-based tests on random inputs:

- **dd = 0** — applying the exterior derivative twice yields zero (cohomological identity)
- **∇·B = 0** — discrete magnetic divergence is identically zero, enforced by construction
- **Circulation conservation** — total circulation over macroscopic loops is preserved for incompressible Euler flows

An operator is not implemented until a property-based test for its invariant exists and passes.

### Design principles

- `comptime` type safety: incompatible operator compositions are compile errors, not runtime failures
- Usability is a correctness property: the type system makes correct code easy, not just incorrect code impossible
- Struct-of-Arrays layout (`std.MultiArrayList`) for all mesh entities
- Explicit allocators everywhere — no hidden heap allocations
- TigerStyle: assert liberally, bound all loops, `const` by default

---

## Development

Development is organized into **epochs** (∼1 month), **milestones** (∼1 week), and **issues** (∼2–6 hours). Each milestone has an explicit mathematical acceptance criterion — it is not done until that invariant passes in CI.

### Workflow

```
/ideate     Pressure-test raw ideas against the vision and architecture;
            produces ideation records and horizon entries

/epoch      Plan a new epoch — back-and-forth conversation producing
            project/epoch_N/roadmap.md and a GitHub Project

/milestone  Given an epoch, create a GitHub Milestone with 10–20 issues,
            labels, and acceptance criterion; adds items to the Project board

/tackle     Pick the highest-priority open issue and work it to completion:
            branch → test-first implementation → PR → CI → squash-merge

/review     Math-aware PR review: checks invariants, comptime type safety,
            SoA layout, TigerStyle, horizons, and process before merge

/status     Live snapshot of current milestone progress — open/closed issues,
            acceptance criterion status, and recommended next action

/decide     Log a non-obvious architectural decision to
            project/epoch_N/decision_log.md
```

These skills are available in [Claude Code](https://claude.ai/claude-code). Run any of them from the project root.

### Issue labels

Issues are organized by a five-prefix taxonomy applied automatically on PRs and manually on issues:

| Prefix | Values |
|--------|--------|
| `type/` | `bug` `feature` `invariant` `perf` `docs` `refactor` `test` `ci` |
| `phase/` | `1` `2` `3` `4` `5` |
| `domain/` | `topology` `forms` `operators` `io` `em` `fluid` `build` |
| `priority/` | `high` `medium` `low` |
| `status/` | `blocked` `needs-decision` `good-first-issue` |

### CI

Five jobs run in parallel on every PR:

| Job | Command |
|-----|---------|
| build | `zig build` |
| test | `zig build test` |
| fmt | `zig fmt --check src/ build.zig` |
| docs | `zig build docs` |
| lint | `zig ast-check` on all source files |

All five must pass before merge. Branch protection is enforced on `main`.

API documentation is deployed to GitHub Pages on every push to `main`.

---

## Contributing

This project is in early development and not yet open to outside contributions. Watch the [Projects](https://github.com/harnesslabs/flux/projects) board for progress.
