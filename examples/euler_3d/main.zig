const euler = @import("euler.zig");

pub fn main() !void {
    return euler.runCli();
}

test {
    _ = @import("euler.zig");
}
