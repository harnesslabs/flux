const std = @import("std");
const runner = @import("runner.zig");

pub fn runExactEvolutionLoop(
    comptime MeshType: type,
    allocator: std.mem.Allocator,
    mesh: *const MeshType,
    evolution: anytype,
    config: runner.RunLoopConfig,
) !runner.RunLoopResult {
    const EvolutionType = @TypeOf(evolution.*);

    const ExactInitializer = struct {
        evolution_ptr: *EvolutionType,

        pub fn fill(self: @This(), _: *const MeshType, values: []f64, time: f64) void {
            self.evolution_ptr.fillExact(time);
            const exact_values = self.evolution_ptr.exactValues();
            std.debug.assert(values.len == exact_values.len);
            if (values.ptr == exact_values.ptr) return;
            @memcpy(values, exact_values);
        }
    };

    const Stepper = struct {
        evolution_ptr: *EvolutionType,

        pub fn step(self: @This(), inner_allocator: std.mem.Allocator) !void {
            try self.evolution_ptr.step(inner_allocator);
        }
    };

    return runner.runExactFieldLoop(
        MeshType,
        allocator,
        mesh,
        evolution.stateValues(),
        evolution.exactValues(),
        config,
        ExactInitializer{ .evolution_ptr = evolution },
        Stepper{ .evolution_ptr = evolution },
    );
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
