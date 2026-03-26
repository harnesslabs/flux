Audit the codebase for style and hygiene issues. Run on `main` between epochs or mid-epoch.

Usage: /audit-style [component]

If a component is specified, scope to that component per `project/components.md`. Otherwise, audit the full `src/` tree.

## What to look for

### Naming
- Full words, no abbreviations: `exterior_derivative` not `ext_deriv`
- Units and qualifiers appended in descending significance
- Consistent naming patterns across similar constructs

### Dead code and stale references
- Unreachable code paths
- Commented-out code blocks
- Inline comments referencing issue numbers, PR numbers, or other artifacts that may be stale
- TODO/FIXME/HACK comments that should be filed as issues or resolved
- Unused imports, variables, or functions

### Comment quality
- Comments should explain *why*, not *what*
- Mathematical content should use Unicode symbols (∂, ∇, ∈, ★, ∧, ⟨⟩, Ω)
- Citations should include author, title, and equation number
- No fabricated URLs

### Code duplication
- Similar logic repeated in multiple places that could share a common implementation
- Copy-paste patterns with minor variations

### Structural consistency
- Similar constructs should follow similar patterns (e.g., all operator types should have the same method shape)
- Public API should be consistent in style (all snake_case, consistent parameter ordering)

### Formatting
- `zig fmt` compliance (should be caught by CI, but verify)
- Consistent file organization (imports, types, functions, tests)

## Output format

```
## Style Audit: <scope>

Date: YYYY-MM-DD

### Stale references
- <file:line> — <what's stale and why>

### Dead code
- <file:line> — <what's dead>

### Naming issues
- <file:line> — <current name> → <suggested name>

### Duplication
- <file1:line> and <file2:line> — <description of duplicated logic>

### Comment issues
- <file:line> — <what's wrong with the comment>

### Summary
<overall assessment, recommended issues to create>
```

For each finding, propose whether it should be:
- **Auto-fixed now** (rename, remove dead code, fix a comment — no design choice)
- **Filed as an issue** (duplication that needs a shared abstraction, structural inconsistency)

Ask the user which findings to act on before making changes or creating issues.
