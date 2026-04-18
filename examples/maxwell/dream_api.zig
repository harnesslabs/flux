//! Dream Maxwell API sketch.
//!
//! This file is intentionally aspirational and is not wired into the build.
//! The point is to pressure-test the public nouns and verbs we actually want.
//!
//! The main target shape is:
//! - `System` is the owned runtime noun.
//! - `Evolution` owns execution policy and attached listeners.
//! - generic methods such as `Leapfrog` validate a comptime contract on the
//!   chosen system instead of re-owning problem state.
//! - fields are structural; diagnostics and measurements are opt-in system
//!   capabilities.

const std = @import("std");
const flux = @import("flux");

const Geometry = enum {
    cartesian,
    spherical,
    toroidal,
};

pub fn runDipole(
    comptime dim: u8,
    comptime geometry: Geometry,
    allocator: std.mem.Allocator,
    writer: anytype,
) !void {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Maxwell = flux.systems.Maxwell(dim, Mesh);

    var mesh = try makeMesh(dim, geometry, allocator, .{
        .counts = countsFor(dim, 64),
        .extents = extentsFor(dim, 1.0),
    });
    defer mesh.deinit(allocator);

    var system = try Maxwell.dipole(allocator, &mesh, .{
        .center = centerFor(dim, 0.5),
        .frequency_hz = 0.5,
        .amplitude = 1.0,
        .boundary = .pec,
    });
    defer system.deinit(allocator);

    var evolution = flux.evolution.Evolution(Maxwell, flux.integrators.Leapfrog)
        .config()
        .dt(system.stableDt())
        .steps(1000)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(Maxwell)
                .field(.electric)
                .field(.magnetic)
                .directory("output/maxwell")
                .baseName("dipole")
                .everySteps(10),
        )
        .listen(
            flux.listeners.Snapshots(Maxwell)
                .measurement(.energy_density)
                .directory("output/maxwell")
                .baseName("dipole_energy_density")
                .everySteps(10),
        )
        .init(allocator, &system);
    defer evolution.deinit();

    try evolution.run();
}

pub fn runCavity(
    comptime dim: u8,
    comptime geometry: Geometry,
    allocator: std.mem.Allocator,
    writer: anytype,
    options: CavityOptions(dim),
) !void {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Maxwell = flux.systems.Maxwell(dim, Mesh);

    var mesh = try makeMesh(dim, geometry, allocator, .{
        .counts = countsFor(dim, 48),
        .extents = extentsFor(dim, 1.0),
    });
    defer mesh.deinit(allocator);

    var system = try Maxwell.cavity(allocator, &mesh, .{
        .mode = options.mode,
        .boundary = options.boundary,
    });
    defer system.deinit(allocator);

    var evolution = flux.evolution.Evolution(Maxwell, flux.integrators.Leapfrog)
        .config()
        .dt(system.stableDt())
        .steps(2000)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(Maxwell)
                .field(.electric)
                .field(.magnetic)
                .measurement(.energy)
                .directory("output/maxwell")
                .baseName("cavity")
                .everySteps(25),
        )
        .init(allocator, &system);
    defer evolution.deinit();

    try evolution.run();
}

pub fn runCavityWithReference(
    comptime dim: u8,
    comptime geometry: Geometry,
    allocator: std.mem.Allocator,
    writer: anytype,
    options: CavityOptions(dim),
) !void {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Maxwell = flux.systems.Maxwell(dim, Mesh);

    var mesh = try makeMesh(dim, geometry, allocator, .{
        .counts = countsFor(dim, 32),
        .extents = extentsFor(dim, 1.0),
    });
    defer mesh.deinit(allocator);

    var system = try Maxwell.cavity(allocator, &mesh, .{
        .mode = options.mode,
        .boundary = options.boundary,
    });
    defer system.deinit(allocator);

    var evolution = flux.evolution.Evolution(Maxwell, flux.integrators.Leapfrog)
        .config()
        .dt(system.stableDt())
        .steps(500)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(Maxwell)
                .field(.electric)
                .field(.magnetic)
                .measurement(.electric_l2)
                .measurement(.magnetic_l2)
                .directory("output/maxwell")
                .baseName("cavity_error")
                .everySteps(10),
        )
        .init(allocator, &system);
    defer evolution.deinit();

    var reference = flux.analysis.ErrorProbe(Maxwell)
        .against(flux.reference.MaxwellCavityMode(dim).init(.{
            .mode = options.mode,
            .extents = extentsFor(dim, 1.0),
        }))
        .measurement(.electric_l2)
        .measurement(.magnetic_l2);

    try evolution.runWith(reference);
}

pub fn runAdaptiveDipole(
    comptime dim: u8,
    comptime geometry: Geometry,
    allocator: std.mem.Allocator,
    writer: anytype,
    options: AdaptiveDipoleOptions(dim),
) !void {
    const Mesh = flux.topology.Mesh(dim, dim);
    const Maxwell = flux.systems.Maxwell(dim, Mesh);
    const Adaptive = flux.integrators.PIDControlledRK(Maxwell);

    var mesh = try makeMesh(dim, geometry, allocator, .{
        .counts = countsFor(dim, options.resolution),
        .extents = extentsFor(dim, options.domain_extent),
    });
    defer mesh.deinit(allocator);

    var system = try Maxwell.dipole(allocator, &mesh, .{
        .center = centerFor(dim, options.source_center),
        .frequency_hz = options.frequency_hz,
        .amplitude = options.amplitude,
        .boundary = options.boundary,
    });
    defer system.deinit(allocator);

    var probe = flux.analysis.ErrorProbe(Maxwell)
        .against(flux.reference.MaxwellDipole(dim).init(.{
            .center = centerFor(dim, options.source_center),
            .frequency_hz = options.frequency_hz,
            .amplitude = options.amplitude,
            .extents = extentsFor(dim, options.domain_extent),
        }))
        .measurement(.electric_l2);

    var evolution = flux.evolution.Evolution(Maxwell, Adaptive)
        .config()
        .dt(options.dt_initial)
        .steps(options.steps_max)
        .listen(flux.listeners.Progress(writer))
        .listen(
            flux.listeners.Snapshots(Maxwell)
                .field(.electric)
                .measurement(.electric_l2)
                .measurement(.time_step)
                .directory("output/maxwell")
                .baseName("adaptive_dipole")
                .everySteps(10),
        )
        .methodOptions(.{
            .tolerance = options.tolerance,
            .dt_min = options.dt_min,
            .dt_max = options.dt_max,
            .safety_factor = 0.9,
        })
        .init(allocator, &system);
    defer evolution.deinit();

    try evolution.runWith(probe);
}

fn makeMesh(
    comptime dim: u8,
    comptime geometry: Geometry,
    allocator: std.mem.Allocator,
    options: anytype,
) !flux.topology.Mesh(dim, dim) {
    const Mesh = flux.topology.Mesh(dim, dim);
    return switch (geometry) {
        .cartesian => Mesh.uniform(allocator, options.counts, options.extents),
        .spherical => Mesh.spherical_shell(allocator, options.counts, .{
            .inner_radius = 0.25,
            .outer_radius = 1.0,
        }),
        .toroidal => Mesh.toroidal_shell(allocator, options.counts, .{
            .major_radius = 1.0,
            .minor_radius = 0.2,
        }),
    };
}

fn countsFor(comptime dim: u8, value: u32) std.meta.Vector(dim, u32) {
    return @splat(value);
}

fn extentsFor(comptime dim: u8, value: f64) std.meta.Vector(dim, f64) {
    return @splat(value);
}

fn centerFor(comptime dim: u8, value: f64) std.meta.Vector(dim, f64) {
    return @splat(value);
}

fn CavityOptions(comptime dim: u8) type {
    return struct {
        mode: flux.reference.Mode(dim),
        boundary: BoundaryCondition = .pec,
    };
}

fn AdaptiveDipoleOptions(comptime dim: u8) type {
    return struct {
        resolution: u32 = 64,
        domain_extent: f64 = 1.0,
        source_center: f64 = 0.5,
        frequency_hz: f64 = 0.5,
        amplitude: f64 = 1.0,
        boundary: BoundaryCondition = .pec,
        dt_initial: f64 = 1.0e-3,
        dt_min: f64 = 1.0e-5,
        dt_max: f64 = 1.0e-2,
        tolerance: f64 = 1.0e-6,
        steps_max: u32 = 5000,
    };
    _ = dim;
}

const BoundaryCondition = enum {
    pec,
    pmc,
};

// Dream system contract sketch:
//
// pub fn Maxwell(comptime dim: u8, comptime MeshType: type) type {
//     return struct {
//         pub const Field = enum { electric, magnetic, current };
//         pub const Diagnostic = enum { energy };
//         pub const Measurement = enum {
//             energy,
//             energy_density,
//             electric_l2,
//             magnetic_l2,
//             time_step,
//         };
//
//         pub fn dipole(allocator: std.mem.Allocator, mesh: *const MeshType, options: DipoleOptions(dim)) !@This() { ... }
//         pub fn cavity(allocator: std.mem.Allocator, mesh: *const MeshType, options: CavityOptions(dim)) !@This() { ... }
//         pub fn stableDt(self: *@This()) f64 { ... }
//         pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void { ... }
//
//         // Used by listeners and output when a field is requested.
//         pub fn field(self: *const @This(), which: Field) FieldView { ... }
//
//         // Available only if the system declares Diagnostics.
//         pub fn diagnostic(self: *const @This(), allocator: std.mem.Allocator, which: Diagnostic) !DiagnosticValue { ... }
//
//         // Available when the system and any attached comparison machinery can
//         // provide the requested quantity.
//         pub fn measurement(self: *const @This(), allocator: std.mem.Allocator, which: Measurement) !MeasurementValue { ... }
//
//         // Optional contract consumed by splitting methods such as Leapfrog.
//         pub const Splitting = struct {
//             pub fn first(self: *@This(), allocator: std.mem.Allocator, dt: f64) !void { ... }
//             pub fn second(self: *@This(), allocator: std.mem.Allocator, dt: f64) !void { ... }
//             pub fn applyBoundary(self: *@This()) void { ... }
//         };
//     };
// }
//
// Then generic methods become capability-driven:
//
// pub fn Leapfrog(comptime System: type) type {
//     comptime flux.concepts.SplitSystem(System);
//     return struct {
//         pub const Options = struct {};
//
//         pub fn advance(allocator: std.mem.Allocator, system: *System, dt: f64, _: Options) !void {
//             try system.splitting.first(allocator, dt);
//             try system.splitting.second(allocator, dt);
//             if (@hasDecl(@TypeOf(system.splitting), "applyBoundary")) {
//                 system.splitting.applyBoundary();
//             }
//         }
//     };
// }

// Desired API properties:
//
// 1. Dimension is mostly an adjective. The same run functions work for 2D/3D.
// 2. Mesh branching is about geometry families, not raw dimension plumbing.
// 3. `System` is the real runtime noun: fields, mesh, operators, boundary data,
//    sources, and caches all live there.
// 4. `Evolution` owns execution policy only: dt, time horizon, listeners.
// 5. `Method` stays thin. The leapfrog noun should not re-own system state.
// 6. Systems expose fields by default and may opt into named diagnostics and
//    measurements through a single capability story.
// 7. Listeners snapshot fields and measurements through one common
//    event surface instead of inventing family-local renderer structs.
// 8. Reference/error machinery is light. `ErrorProbe` composes with evolution
//    and should also be reusable for AMR, adaptivity, and steady-state studies.
// 9. Adaptive methods still fit this model: method-specific runtime control
//    state belongs to the method/options layer, not to the base `Evolution`.
