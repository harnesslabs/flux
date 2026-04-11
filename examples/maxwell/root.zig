const impl = @import("impl.zig");

pub const Mesh = impl.Mesh;
pub const Config = impl.Config;
pub const RunResult = impl.RunResult;
pub const State = impl.State;
pub const run = impl.run;
pub const makeMesh = impl.makeMesh;
pub const step = impl.step;
pub const seedReferenceMode = impl.seedReferenceMode;
pub const divergenceNorm = impl.divergenceNorm;

test {
    _ = impl;
    _ = @import("tests.zig");
}
