const runner = @import("runner.zig");

pub const RunLoopConfig = runner.RunLoopConfig;
pub const RunLoopResult = runner.RunLoopResult;
pub const ExactFieldRenderer = runner.ExactFieldRenderer;
pub const runConvergenceStudy = runner.runConvergenceStudy;
pub const runExactFieldLoop = runner.runExactFieldLoop;
pub const runSimulationLoop = runner.runSimulationLoop;

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
