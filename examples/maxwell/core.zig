const std = @import("std");
const flux = @import("flux");

const cochain = flux.forms;
const operator_context_mod = flux.operators.context;

pub fn Core(comptime MeshType: type, comptime Hooks: type) type {
    return struct {
        pub const State = struct {
            const Self = @This();

            pub const OneForm = cochain.Cochain(MeshType, 1, cochain.Primal);
            pub const TwoForm = cochain.Cochain(MeshType, 2, cochain.Primal);

            E: OneForm,
            B: TwoForm,
            J: OneForm,
            mesh: *const MeshType,
            operators: *operator_context_mod.OperatorContext(MeshType),
            timestep: u64,

            pub fn init(allocator: std.mem.Allocator, mesh: *const MeshType) !Self {
                var electric = try OneForm.init(allocator, mesh);
                errdefer electric.deinit(allocator);

                var magnetic = try TwoForm.init(allocator, mesh);
                errdefer magnetic.deinit(allocator);

                var current = try OneForm.init(allocator, mesh);
                errdefer current.deinit(allocator);

                const operators = try operator_context_mod.OperatorContext(MeshType).init(allocator, mesh);
                errdefer operators.deinit();
                try Hooks.primeOperators(operators);

                return .{
                    .E = electric,
                    .B = magnetic,
                    .J = current,
                    .mesh = mesh,
                    .operators = operators,
                    .timestep = 0,
                };
            }

            pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
                self.operators.deinit();
                self.J.deinit(allocator);
                self.B.deinit(allocator);
                self.E.deinit(allocator);
            }
        };

        pub fn faradayStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            var derivative = try (try state.operators.exteriorDerivative(cochain.Primal, 1)).apply(allocator, state.E);
            defer derivative.deinit(allocator);

            derivative.scale(dt);
            state.B.sub(derivative);
        }

        pub fn applyPecBoundary(state: anytype) void {
            for (state.mesh.boundary_edges) |edge_idx| {
                state.E.values[edge_idx] = 0.0;
            }
        }

        pub fn ampereStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            try Hooks.ampereStep(allocator, state, dt);
        }

        pub fn leapfrogStep(allocator: std.mem.Allocator, state: anytype, dt: f64) !void {
            try faradayStep(allocator, state, dt);
            try ampereStep(allocator, state, dt);
            if (Hooks.apply_boundary_in_leapfrog) {
                applyPecBoundary(state);
            }
            state.timestep += 1;
        }
    };
}
