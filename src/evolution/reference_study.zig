//! Reference-study helpers for pairing an evolving system with an analytic or
//! otherwise known comparison field.

const std = @import("std");
const testing = std.testing;

pub fn ReferenceStudy(
    comptime MeshType: type,
    comptime InitializerType: type,
    comptime ErrorMeasureType: type,
) type {
    return struct {
        mesh: *const MeshType,
        state_values: []const f64,
        exact_values: []f64,
        initializer: InitializerType,
        error_measure: ErrorMeasureType,

        pub fn init(
            allocator: std.mem.Allocator,
            mesh: *const MeshType,
            state_values: []const f64,
            initializer: InitializerType,
            error_measure: ErrorMeasureType,
        ) !@This() {
            return .{
                .mesh = mesh,
                .state_values = state_values,
                .exact_values = try allocator.alloc(f64, state_values.len),
                .initializer = initializer,
                .error_measure = error_measure,
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            allocator.free(self.exact_values);
        }

        pub fn fillExact(self: *@This(), time: f64) void {
            self.initializer.fill(self.mesh, self.exact_values, time);
        }

        pub fn exactValues(self: *@This()) []const f64 {
            return self.exact_values;
        }

        pub fn stateValues(self: *const @This()) []const f64 {
            return self.state_values;
        }

        pub fn l2Error(self: *const @This()) f64 {
            return self.error_measure.compute(self.mesh, self.state_values, self.exact_values);
        }
    };
}

pub fn empiricalRates(
    allocator: std.mem.Allocator,
    errors: []const f64,
    refinement_ratio: f64,
) ![]f64 {
    std.debug.assert(refinement_ratio > 1.0);

    if (errors.len < 2) {
        return allocator.alloc(f64, 0);
    }

    const rates = try allocator.alloc(f64, errors.len - 1);
    errdefer allocator.free(rates);

    for (0..errors.len - 1) |idx| {
        const coarse = errors[idx];
        const fine = errors[idx + 1];
        std.debug.assert(coarse > 0.0);
        std.debug.assert(fine > 0.0);
        rates[idx] = std.math.log(f64, refinement_ratio, coarse / fine);
    }

    return rates;
}

const MockMesh = struct {};

const MockExactInitializer = struct {
    pub fn fill(_: @This(), _: *const MockMesh, values: []f64, time: f64) void {
        for (values, 0..) |*value, idx| {
            value.* = time + @as(f64, @floatFromInt(idx));
        }
    }
};

const MockErrorMeasure = struct {
    pub fn compute(_: @This(), _: *const MockMesh, approx: []const f64, exact: []const f64) f64 {
        var sum: f64 = 0.0;
        for (approx, exact) |approx_value, exact_value| {
            sum += @abs(approx_value - exact_value);
        }
        return sum;
    }
};

test "ReferenceStudy owns exact buffer and computes error generically" {
    const allocator = testing.allocator;
    var mesh = MockMesh{};
    var state = [_]f64{ 1.0, 2.0 };

    const Study = ReferenceStudy(MockMesh, MockExactInitializer, MockErrorMeasure);
    var study = try Study.init(
        allocator,
        &mesh,
        state[0..],
        MockExactInitializer{},
        MockErrorMeasure{},
    );
    defer study.deinit(allocator);

    study.fillExact(0.0);
    try testing.expectEqual(state.len, study.exactValues().len);
    try testing.expectApproxEqAbs(2.0, study.l2Error(), 1e-15);
}

test "empiricalRates returns pairwise convergence rates" {
    const allocator = testing.allocator;
    const errors = [_]f64{ 4.0, 1.0, 0.25 };
    const rates = try empiricalRates(allocator, &errors, 2.0);
    defer allocator.free(rates);

    try testing.expectEqual(@as(usize, 2), rates.len);
    try testing.expectApproxEqAbs(2.0, rates[0], 1e-12);
    try testing.expectApproxEqAbs(2.0, rates[1], 1e-12);
}
