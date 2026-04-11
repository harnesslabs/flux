const std = @import("std");
const runner = @import("runner.zig");

pub fn runEvolutionLoop(
    allocator: std.mem.Allocator,
    evolution: anytype,
    config: runner.RunLoopConfig,
    renderer: anytype,
) !runner.RunLoopResult {
    const Capturer = struct {
        renderer_value: @TypeOf(renderer),

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            try series.capture(time, self.renderer_value);
        }
    };

    return runner.runStepLoop(
        allocator,
        config,
        evolution,
        Capturer{ .renderer_value = renderer },
    );
}

pub fn runExactEvolutionLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    evolution: anytype,
    config: runner.RunLoopConfig,
) !runner.RunLoopResult {
    const EvolutionType = @TypeOf(evolution.*);

    const Capturer = struct {
        mesh_ptr: *const MeshType,
        evolution_ptr: *EvolutionType,

        pub fn capture(self: @This(), series: anytype, time: f64) !void {
            self.evolution_ptr.fillExact(time);
            try series.capture(time, runner.ExactFieldRenderer(MeshType){
                .mesh = self.mesh_ptr,
                .state = self.evolution_ptr.stateValues(),
                .exact = self.evolution_ptr.exactValues(),
            });
        }
    };

    const result = try runner.runStepLoop(
        allocator,
        config,
        evolution,
        Capturer{
            .mesh_ptr = mesh,
            .evolution_ptr = evolution,
        },
    );
    evolution.fillExact(config.final_time);
    return result;
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
