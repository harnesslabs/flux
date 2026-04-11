const impl = @import("impl.zig");

pub const Demo = impl.Demo;
pub const Config2D = impl.Config2D;
pub const RunResult2D = impl.RunResult2D;
pub const Config3D = impl.Config3D;
pub const RunResult3D = impl.RunResult3D;
pub const MaxwellState3D = impl.MaxwellState3D;

pub const StateForMesh2D = impl.StateForMesh2D;
pub const StateForMesh3D = impl.StateForMesh3D;
pub const Mesh = impl.Mesh;
pub const Config = impl.Config;
pub const RunResult = impl.RunResult;
pub const State = impl.State;

pub const MaxwellSystem = impl.MaxwellSystem;
pub const MaxwellLeapfrog = impl.MaxwellLeapfrog;
pub const PointDipole = impl.PointDipole;

pub const faraday_step = impl.faraday_step;
pub const ampere_step = impl.ampere_step;
pub const apply_pec_boundary = impl.apply_pec_boundary;
pub const leapfrog_step = impl.leapfrog_step;
pub const electromagnetic_energy = impl.electromagnetic_energy;

pub const project_te10_e = impl.project_te10_e;
pub const project_te10_b = impl.project_te10_b;
pub const faradayStep3D = impl.faradayStep3D;
pub const ampereStep3D = impl.ampereStep3D;
pub const applyPecBoundary3D = impl.applyPecBoundary3D;
pub const leapfrogStep3D = impl.leapfrogStep3D;
pub const runSimulation = impl.runSimulation;
pub const project_tm110_e = impl.project_tm110_e;
pub const project_tm110_potential = impl.project_tm110_potential;
pub const project_tm110_b = impl.project_tm110_b;
pub const seedTm110Mode = impl.seedTm110Mode;
pub const writeSnapshot = impl.writeSnapshot;
pub const makeCavityMesh = impl.makeCavityMesh;
pub const divergenceNorm3D = impl.divergenceNorm3D;

pub const run = impl.run;
pub const runDriver = impl.runDriver;
pub const makeMesh = impl.makeMesh;
pub const leapfrogStep = impl.leapfrogStep;
pub const step = impl.step;
pub const seedReferenceMode = impl.seedReferenceMode;
pub const divergenceNorm = impl.divergenceNorm;
pub const projectReferenceMagneticField = impl.projectReferenceMagneticField;

test {
    _ = impl;
}
