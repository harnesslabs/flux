const impl = @import("impl.zig");

pub const Mesh2D = impl.Mesh2D;
pub const VertexVorticity = impl.VertexVorticity;
pub const FaceVorticity = impl.FaceVorticity;
pub const FaceTracer = impl.FaceTracer;
pub const Demo = impl.Demo;
pub const Config2D = impl.Config2D;
pub const RunResult2D = impl.RunResult2D;
pub const State2D = impl.State2D;

pub const Mesh3D = impl.Mesh3D;
pub const Velocity = impl.Velocity;
pub const Vorticity = impl.Vorticity;
pub const Potential = impl.Potential;
pub const Config3D = impl.Config3D;
pub const RunResult3D = impl.RunResult3D;
pub const State3D = impl.State3D;

pub const initializeGaussianVortex = impl.initializeGaussianVortex;
pub const initializeVortexDipole = impl.initializeVortexDipole;
pub const totalCirculation = impl.totalCirculation;
pub const step2D = impl.step2D;
pub const run2D = impl.run2D;
pub const seedReferenceMode3D = impl.seedReferenceMode3D;
pub const recoverVelocityFromVorticity = impl.recoverVelocityFromVorticity;
pub const step3D = impl.step3D;
pub const run3D = impl.run3D;
pub const computeHelicity3D = impl.computeHelicity3D;

pub const Mesh = impl.Mesh;
pub const Config = impl.Config;
pub const RunResult = impl.RunResult;
pub const State = impl.State;
pub const run = impl.run;
pub const step = impl.step;
pub const seedReferenceMode = impl.seedReferenceMode;
pub const conservedQuantity = impl.conservedQuantity;

test {
    _ = impl;
}
