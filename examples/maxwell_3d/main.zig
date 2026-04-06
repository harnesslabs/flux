const maxwell = @import("maxwell.zig");

pub fn main() !void {
    return maxwell.runCli();
}

test {
    _ = @import("maxwell.zig");
}
