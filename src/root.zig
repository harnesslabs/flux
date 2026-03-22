pub const io = @import("io/vtk.zig");
pub const math = @import("math/sparse.zig");
pub const topology = @import("topology/mesh.zig");

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
