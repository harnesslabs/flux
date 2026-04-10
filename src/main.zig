const app = @import("example_app");

pub fn main() !void {
    try app.main();
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
