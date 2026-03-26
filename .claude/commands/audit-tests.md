Audit the test suite for coverage, quality, and redundancy. Run on `main` between epochs or mid-epoch.

Usage: /audit-tests [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

## What to look for

### Coverage gaps
- Public functions without any test coverage
- Operators without property-based tests on random inputs
- Edge cases not tested: empty mesh, single element, boundary-only elements, degenerate geometry
- Error paths: what happens with invalid inputs? Are error returns tested?

### Property test quality
- Property tests should use random inputs, not hand-crafted examples only
- Random inputs should cover a meaningful range (different mesh sizes, different cochain values)
- The number of random samples should be sufficient (≥100 for quick tests, ≥1000 for critical invariants)
- Property tests must test the *mathematical* invariant, not just "it doesn't crash"

### Redundant tests
- Tests that verify the same property in essentially the same way
- Tests that are subsumed by a more general property test
- Tests that test internal implementation details rather than observable behavior

### Test naming
- Names should state the invariant: `test "dd = 0 for random 1-forms on triangular mesh"`
- Names should not be generic: `test "basic test"` or `test "it works"` are not acceptable

### Convergence tests
- For operators with known analytical solutions: does error decrease at the expected rate when the mesh is refined?
- Are convergence rates checked quantitatively (assert rate ≥ expected order - tolerance)?

### Comptime tests
- Are there `comptime` blocks that verify type-level properties?
- Could any runtime test be promoted to a compile-time guarantee?

### Test organization
- Tests should be in the same file as the code they test (Zig convention)
- Test helpers should be clearly marked and reusable
- Shared test fixtures (meshes, cochains) should be factored out if used across multiple test blocks

## Output format

```
## Test Audit: <scope>

Date: YYYY-MM-DD

### Coverage gaps
- <function/module> — <what's not tested, what invariant should be verified>

### Weak tests (exist but insufficient)
- <test name in file> — <why it's weak, what would make it strong>

### Redundant tests
- <test name> and <test name> — <why one subsumes the other>

### Missing convergence tests
- <operator/capability> — <what analytical solution to compare against>

### Comptime promotion opportunities
- <runtime test> — <could be a compile-time guarantee because ...>

### Summary
<N gaps, M weak tests, K redundant, L convergence opportunities>
<overall test suite health assessment>
```

For each finding, propose whether it should be:
- **Auto-fixed now** (rename a test, remove a clearly redundant test)
- **Filed as an issue** (new property test, convergence test, coverage expansion)

Convergence tests and property tests should always be filed as issues — they require careful mathematical thought, not quick fixes.

Ask the user which findings to act on.
