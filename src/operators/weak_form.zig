//! Local-to-global assembly for FEEC weak-form operators.

const std = @import("std");
const sparse = @import("../math/sparse.zig");

pub const AssembleError = error{
    NotYetImplemented,
};

pub fn assemble(
    comptime k: comptime_int,
    allocator: std.mem.Allocator,
    mesh: anytype,
    kernel: anytype,
) !sparse.CsrMatrix(f64) {
    _ = k;
    _ = allocator;
    _ = mesh;
    _ = kernel;
    return AssembleError.NotYetImplemented;
}
