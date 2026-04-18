const std = @import("std");

pub fn Progress(writer: *std.Io.Writer) ProgressListener {
    return .{ .writer = writer };
}

pub const ProgressListener = struct {
    writer: *std.Io.Writer,
    timer: ?std.time.Timer = null,
    last_draw_ns: u64 = 0,
    bar_width: u32 = 40,

    pub fn onRunBegin(self: *@This(), _: anytype) !void {
        self.timer = std.time.Timer.start() catch
            @panic("OS timer unavailable — cannot run simulation");
        self.last_draw_ns = 0;
    }

    pub fn onStep(self: *@This(), event: anytype) !void {
        const timer = if (self.timer) |*timer_value| timer_value else return;
        const elapsed_ns = timer.read();

        if (elapsed_ns - self.last_draw_ns < 50_000_000 and event.step_index < event.step_goal) return;
        self.last_draw_ns = elapsed_ns;

        const elapsed_s = @as(f64, @floatFromInt(elapsed_ns)) / 1_000_000_000.0;
        const frac = @as(f64, @floatFromInt(event.step_index)) / @as(f64, @floatFromInt(event.step_goal));
        const pct = frac * 100.0;
        const steps_per_sec = if (elapsed_s > 0.01) @as(f64, @floatFromInt(event.step_index)) / elapsed_s else 0.0;
        const remaining = @as(f64, @floatFromInt(event.step_goal - event.step_index));
        const eta_s = if (steps_per_sec > 0.01) remaining / steps_per_sec else 0.0;

        const filled: u32 = @intFromFloat(frac * @as(f64, @floatFromInt(self.bar_width)));

        var bar: [64]u8 = undefined;
        for (0..self.bar_width) |idx| {
            bar[idx] = if (idx < filled) '#' else '-';
        }
        const bar_str = bar[0..self.bar_width];

        var elapsed_buf: [16]u8 = undefined;
        var eta_buf: [16]u8 = undefined;
        try self.writer.print("\r  {s}  {d:>5.1}%  {d}/{d}  {s}  ETA {s}  {d:.0} steps/s    ", .{
            bar_str,
            pct,
            event.step_index,
            event.step_goal,
            formatDuration(&elapsed_buf, elapsed_s),
            formatDuration(&eta_buf, eta_s),
            steps_per_sec,
        });
        try self.writer.flush();
    }

    pub fn onRunEnd(self: *@This(), _: anytype) !void {
        try self.writer.writeAll("\r");
        for (0..120) |_| try self.writer.writeByte(' ');
        try self.writer.writeAll("\r");
        try self.writer.flush();
    }
};

pub fn Snapshots(comptime SystemType: type) SnapshotListener(SystemType) {
    return .{};
}

pub fn SnapshotListener(comptime SystemType: type) type {
    const SystemDeclType = switch (@typeInfo(SystemType)) {
        .pointer => |pointer| pointer.child,
        else => SystemType,
    };

    if (!@hasDecl(SystemDeclType, "Field")) {
        @compileError("Snapshots requires the system to declare `pub const Field`");
    }

    const Field = SystemDeclType.Field;
    const Measurement = if (@hasDecl(SystemDeclType, "Measurement")) SystemDeclType.Measurement else void;
    const MeasurementStorage = if (Measurement == void) u8 else Measurement;
    const field_capacity = 8;
    const measurement_capacity = 8;

    return struct {
        const Self = @This();

        fields: [field_capacity]Field = undefined,
        field_count: usize = 0,
        measurements: [measurement_capacity]MeasurementStorage = undefined,
        measurement_count: usize = 0,
        output_dir: []const u8 = "",
        base_name: []const u8 = "snapshot",
        interval_steps: u32 = 1,
        capture_initial: bool = true,
        emit_final: bool = true,
        last_captured_step: ?u32 = null,
        csv_initialized: bool = false,

        pub fn field(self: Self, which: Field) Self {
            std.debug.assert(self.field_count < field_capacity);
            var next = self;
            next.fields[next.field_count] = which;
            next.field_count += 1;
            return next;
        }

        pub fn measurement(self: Self, which: Measurement) Self {
            if (Measurement == void) {
                @compileError("Snapshots.measurement requires the chosen system to declare `pub const Measurement`");
            }
            std.debug.assert(self.measurement_count < measurement_capacity);
            var next = self;
            next.measurements[next.measurement_count] = which;
            next.measurement_count += 1;
            return next;
        }

        pub fn directory(self: Self, value: []const u8) Self {
            var next = self;
            next.output_dir = value;
            return next;
        }

        pub fn baseName(self: Self, value: []const u8) Self {
            var next = self;
            next.base_name = value;
            return next;
        }

        pub fn everySteps(self: Self, value: u32) Self {
            std.debug.assert(value > 0);
            var next = self;
            next.interval_steps = value;
            return next;
        }

        pub fn onRunBegin(self: *Self, event: anytype) !void {
            if (!self.capture_initial) return;
            try self.capture(event);
        }

        pub fn onStep(self: *Self, event: anytype) !void {
            if (event.step_index % self.interval_steps != 0) return;
            try self.capture(event);
        }

        pub fn onRunEnd(self: *Self, event: anytype) !void {
            if (!self.emit_final) return;
            if (self.last_captured_step != null and self.last_captured_step.? == event.step_index) return;
            try self.capture(event);
        }

        fn capture(self: *Self, event: anytype) !void {
            std.debug.assert(self.output_dir.len > 0);
            try std.fs.cwd().makePath(self.output_dir);

            if (self.field_count > 0) {
                try self.captureFields(event);
            }
            if (self.measurement_count > 0) {
                try self.captureMeasurements(event);
            }

            self.last_captured_step = event.step_index;
        }

        fn captureFields(self: *Self, event: anytype) !void {
            var path_buf: [std.fs.max_path_bytes]u8 = undefined;
            const path = try std.fmt.bufPrint(&path_buf, "{s}/{s}_{d:0>5}.vtu", .{
                self.output_dir,
                self.base_name,
                event.step_index,
            });

            const file = try std.fs.cwd().createFile(path, .{ .truncate = true });
            defer file.close();

            var buffer: [4096]u8 = undefined;
            var file_writer = file.writer(&buffer);
            try event.system.writeFields(std.heap.page_allocator, &file_writer.interface, self.fields[0..self.field_count]);
        }

        fn captureMeasurements(self: *Self, event: anytype) !void {
            var path_buf: [std.fs.max_path_bytes]u8 = undefined;
            const path = try std.fmt.bufPrint(&path_buf, "{s}/{s}_measurements.csv", .{
                self.output_dir,
                self.base_name,
            });

            const file = if (!self.csv_initialized)
                try std.fs.cwd().createFile(path, .{ .truncate = true })
            else
                try std.fs.cwd().openFile(path, .{ .mode = .write_only });
            defer file.close();

            if (self.csv_initialized) {
                try file.seekFromEnd(0);
            }

            var buffer: [4096]u8 = undefined;
            var file_writer = file.writer(&buffer);
            const writer = &file_writer.interface;
            if (!self.csv_initialized) {
                try writer.writeAll("step,time");
                for (self.measurements[0..self.measurement_count]) |measurement_value| {
                    const selected_measurement = @as(Measurement, measurement_value);
                    try writer.print(",{s}", .{@tagName(selected_measurement)});
                }
                try writer.writeByte('\n');
                self.csv_initialized = true;
            }

            try writer.print("{d},{d:.12}", .{ event.step_index, event.time });
            for (self.measurements[0..self.measurement_count]) |measurement_value| {
                const selected_measurement = @as(Measurement, measurement_value);
                const value = try event.system.measurement(std.heap.page_allocator, selected_measurement);
                try writer.print(",{d:.16}", .{value});
            }
            try writer.writeByte('\n');
        }
    };
}

fn formatDuration(buf: *[16]u8, seconds: f64) []const u8 {
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
