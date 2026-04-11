const impl = @import("impl.zig");

pub const SurfaceKind = impl.SurfaceKind;

pub const Mesh = impl.Mesh;
pub const State = impl.State;
pub const Config = impl.Config;
pub const RunResult = impl.RunResult;
pub const ConvergenceResult = impl.ConvergenceResult;
pub const run = impl.run;
pub const runConvergenceStudy = impl.runConvergenceStudy;

test {
    _ = impl;
    _ = @import("tests.zig");
}
