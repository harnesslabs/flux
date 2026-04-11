const impl = @import("impl.zig");

pub const SurfaceKind = impl.SurfaceKind;

pub const SurfaceMesh = impl.SurfaceMesh;
pub const SphereVertexField = impl.SphereVertexField;
pub const SphereEdgeField = impl.SphereEdgeField;
pub const SphereConfigImpl = impl.SphereConfigImpl;
pub const SphereRunResultImpl = impl.SphereRunResultImpl;
pub const SphereConvergenceResultImpl = impl.SphereConvergenceResultImpl;

pub const Mesh2D = impl.Mesh2D;
pub const VertexField = impl.VertexField;
pub const PlaneConfigImpl = impl.PlaneConfigImpl;
pub const PlaneRunResultImpl = impl.PlaneRunResultImpl;
pub const PlaneConvergenceResultImpl = impl.PlaneConvergenceResultImpl;

pub const PlaneConfig = impl.PlaneConfig;
pub const PlaneRunResult = impl.PlaneRunResult;
pub const PlaneConvergenceResult = impl.PlaneConvergenceResult;
pub const SphereConfig = impl.SphereConfig;
pub const SphereRunResult = impl.SphereRunResult;
pub const SphereConvergenceResult = impl.SphereConvergenceResult;

pub const runSphereImpl = impl.runSphereImpl;
pub const runSphereConvergenceStudyImpl = impl.runSphereConvergenceStudyImpl;
pub const runPlane = impl.runPlane;
pub const runPlaneConvergenceStudy = impl.runPlaneConvergenceStudy;
pub const runSphere = impl.runSphere;
pub const runSphereConvergenceStudy = impl.runSphereConvergenceStudy;

pub const Mesh = impl.Mesh;
pub const State = impl.State;
pub const Config = impl.Config;
pub const RunResult = impl.RunResult;
pub const ConvergenceResult = impl.ConvergenceResult;
pub const run = impl.run;
pub const runConvergenceStudy = impl.runConvergenceStudy;

test {
    _ = impl;
}
