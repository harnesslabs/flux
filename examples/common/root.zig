//! Shared infrastructure for the flux example suite.
//!
//! Imported as `examples_common` from the umbrella `flux-examples` binary
//! and from each individual example sub-module. Re-exports the CLI parser,
//! snapshot writer, and progress display so callers only need a single
//! `@import("examples_common")`.

pub const cli = @import("cli.zig");
pub const snapshot = @import("snapshot.zig");
pub const progress = @import("progress.zig");

pub const Common = cli.Common;
pub const Parser = cli.Parser;
pub const ParseError = cli.ParseError;
pub const Plan = snapshot.Plan;
pub const Series = snapshot.Series;
pub const ensureDir = snapshot.ensureDir;
pub const Progress = progress.Progress;
pub const formatDuration = progress.formatDuration;

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
