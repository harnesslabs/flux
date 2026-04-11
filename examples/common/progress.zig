//! Console progress bar shared by every example.
//!
//! Uses a concrete `std.Io.Writer` pointer so examples can pass the same live
//! stderr stream through every family and flush redraws immediately.

const std = @import("std");

pub const Progress = struct {
    const Self = @This();

    writer: *std.Io.Writer,
    total: u32,
    timer: std.time.Timer,
    last_draw_ns: u64 = 0,
    bar_width: u32 = 40,

    pub fn init(writer: *std.Io.Writer, total: u32) Self {
        return .{
            .writer = writer,
            .total = total,
            .timer = std.time.Timer.start() catch
                @panic("OS timer unavailable — cannot run simulation"),
        };
    }

    pub fn update(self: *Self, step: u32) void {
        const elapsed_ns = self.timer.read();

        // Cap redraws at ~20 fps so the bar doesn't dominate stderr on
        // fast machines. Always draw the final frame so the bar lands
        // at 100%.
        if (elapsed_ns - self.last_draw_ns < 50_000_000 and step < self.total) return;
        self.last_draw_ns = elapsed_ns;

        const elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0;
        const frac = @as(f64, @floatFromInt(step)) / @as(f64, @floatFromInt(self.total));
        const pct = frac * 100.0;

        // Guard division by zero in the first few milliseconds.
        const steps_per_sec = if (elapsed_s > 0.01) @as(f64, @floatFromInt(step)) / elapsed_s else 0.0;
        const remaining = @as(f64, @floatFromInt(self.total - step));
        const eta_s = if (steps_per_sec > 0.01) remaining / steps_per_sec else 0.0;

        const filled: u32 = @intFromFloat(frac * @as(f64, @floatFromInt(self.bar_width)));

        var bar: [64]u8 = undefined;
        for (0..self.bar_width) |j| {
            bar[j] = if (j < filled) '#' else '-';
        }
        const bar_str = bar[0..self.bar_width];

        var elapsed_buf: [16]u8 = undefined;
        var eta_buf: [16]u8 = undefined;
        self.writer.print("\r  {s}  {d:>5.1}%  {d}/{d}  {s}  ETA {s}  {d:.0} steps/s    ", .{
            bar_str,
            pct,
            step,
            self.total,
            formatDuration(&elapsed_buf, elapsed_s),
            formatDuration(&eta_buf, eta_s),
            steps_per_sec,
        }) catch return;
        self.writer.flush() catch return;
    }

    pub fn finish(self: *Self) void {
        self.writer.writeAll("\r") catch return;
        for (0..120) |_| self.writer.writeByte(' ') catch return;
        self.writer.writeAll("\r") catch return;
        self.writer.flush() catch return;
    }

    pub fn elapsed(self: *Self) f64 {
        return @as(f64, @floatFromInt(self.timer.read())) / 1_000_000_000.0;
    }
};

pub fn formatDuration(buf: *[16]u8, seconds: f64) []const u8 {
    if (seconds < 60.0) {
        return std.fmt.bufPrint(buf, "{d:.1}s", .{seconds}) catch "??";
    }
    const mins: u32 = @intFromFloat(seconds / 60.0);
    const secs: u32 = @intFromFloat(@mod(seconds, 60.0));
    return std.fmt.bufPrint(buf, "{d}m{d:0>2}s", .{ mins, secs }) catch "??";
}

const testing = std.testing;

test "formatDuration sub-minute prints seconds" {
    var buf: [16]u8 = undefined;
    try testing.expectEqualStrings("12.3s", formatDuration(&buf, 12.34));
}

test "formatDuration over a minute prints mins+secs" {
    var buf: [16]u8 = undefined;
    try testing.expectEqualStrings("1m05s", formatDuration(&buf, 65.0));
}

test "formatDuration zero" {
    var buf: [16]u8 = undefined;
    try testing.expectEqualStrings("0.0s", formatDuration(&buf, 0.0));
}
