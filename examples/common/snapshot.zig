//! Time-series snapshot bookkeeping for the flux example suite.
//!
//! Every example simulates for some number of timesteps and writes a small
//! number of VTK frames at uniform intervals, plus a `.pvd` collection that
//! ParaView reads as an animation. The interval math, filename buffering,
//! directory creation, and PVD finalization are identical across examples,
//! so this module owns the pattern and exposes a single `Series` type that
//! examples drive with example-specific render closures.
//!
//! The renderer is supplied by the caller as `anytype`, so this module makes
//! no assumptions about which fields a given simulation writes — only about
//! the bookkeeping around them.

const std = @import("std");
const flux = @import("flux");

const flux_io = flux.io;

/// Plan describing how often to capture snapshots and how many to budget for.
/// `interval == 0` disables snapshot writing entirely (capacity is then 0).
pub const Plan = struct {
    interval: u32,
    capacity: u32,

    /// Compute an evenly-spaced snapshot plan covering `total_steps` with at
    /// most `frames` frames. Returns a disabled plan when either bound is
    /// zero. The capacity is `total_steps / interval + 1` to account for the
    /// optional initial frame; the `+|` saturating add guards the
    /// pathological case `total_steps == u32_max && interval == 1`.
    pub fn fromTotal(total_steps: u32, frames: u32) Plan {
        if (frames == 0 or total_steps == 0) return .{ .interval = 0, .capacity = 0 };
        const interval: u32 = @max(@as(u32, 1), total_steps / frames);
        const capacity: u32 = (total_steps / interval) +| 1;
        return .{ .interval = interval, .capacity = capacity };
    }
};

/// Mutable snapshot collector. Owns the PVD entry buffer and the filename
/// buffer pool, both sized at construction time so the simulation hot loop
/// never allocates.
pub const Series = struct {
    allocator: std.mem.Allocator,
    output_dir: []const u8,
    base_name: []const u8,
    plan: Plan,
    pvd_entries: []flux_io.PvdEntry,
    filename_bufs: [][flux_io.max_snapshot_filename_length]u8,
    count: u32,

    pub fn init(
        allocator: std.mem.Allocator,
        output_dir: []const u8,
        base_name: []const u8,
        plan: Plan,
    ) !Series {
        std.debug.assert(base_name.len > 0);
        std.debug.assert(base_name.len <= 200);
        if (plan.capacity == 0) {
            return .{
                .allocator = allocator,
                .output_dir = output_dir,
                .base_name = base_name,
                .plan = plan,
                .pvd_entries = &.{},
                .filename_bufs = &.{},
                .count = 0,
            };
        }
        try ensureDir(output_dir);
        const pvd_entries = try allocator.alloc(flux_io.PvdEntry, plan.capacity);
        errdefer allocator.free(pvd_entries);
        const filename_bufs = try allocator.alloc(
            [flux_io.max_snapshot_filename_length]u8,
            plan.capacity,
        );
        return .{
            .allocator = allocator,
            .output_dir = output_dir,
            .base_name = base_name,
            .plan = plan,
            .pvd_entries = pvd_entries,
            .filename_bufs = filename_bufs,
            .count = 0,
        };
    }

    pub fn deinit(self: *Series) void {
        if (self.plan.capacity == 0) return;
        self.allocator.free(self.pvd_entries);
        self.allocator.free(self.filename_bufs);
    }

    pub fn enabled(self: Series) bool {
        return self.plan.interval > 0;
    }

    /// Should a snapshot be captured at this 1-based step index?
    /// Returns false when the plan is disabled or the budget is full.
    pub fn dueAt(self: Series, step_one_based: u32) bool {
        if (!self.enabled()) return false;
        if (self.count >= self.plan.capacity) return false;
        return step_one_based % self.plan.interval == 0;
    }

    /// Capture a frame at simulation time `time`. The `renderer` is any
    /// value with a `render(allocator, writer) !void` method. Zig has no
    /// first-class closures, so callers package the data they need to
    /// serialize into a small struct with a render method:
    ///
    ///     const Renderer = struct {
    ///         state: *const State,
    ///         pub fn render(self: @This(), a: Allocator, w: anytype) !void {
    ///             try flux.io.write(w, ..., self.state.mesh.*, ...);
    ///         }
    ///     };
    ///     try series.capture(time, Renderer{ .state = &state });
    ///
    /// The serialized buffer is then written to disk and recorded in the
    /// PVD entry list.
    pub fn capture(
        self: *Series,
        time: f64,
        renderer: anytype,
    ) !void {
        std.debug.assert(self.enabled());
        std.debug.assert(self.count < self.plan.capacity);

        const filename = flux_io.snapshot_filename(
            &self.filename_bufs[self.count],
            self.base_name,
            self.count,
        );

        var output: std.ArrayListUnmanaged(u8) = .{};
        defer output.deinit(self.allocator);
        try renderer.render(self.allocator, output.writer(self.allocator));

        var dir = try std.fs.cwd().openDir(self.output_dir, .{});
        defer dir.close();
        const file = try dir.createFile(filename, .{});
        defer file.close();
        try file.writeAll(output.items);

        self.pvd_entries[self.count] = .{
            .timestep = time,
            .filename = filename,
        };
        self.count += 1;
    }

    /// Write the `.pvd` collection file referencing all captured snapshots.
    /// No-op when nothing was captured.
    pub fn finalize(self: *Series) !void {
        if (self.count == 0) return;

        var output: std.ArrayListUnmanaged(u8) = .{};
        defer output.deinit(self.allocator);
        try flux_io.write_pvd(output.writer(self.allocator), self.pvd_entries[0..self.count]);

        var pvd_buf: [flux_io.max_snapshot_filename_length]u8 = undefined;
        const pvd_name = std.fmt.bufPrint(&pvd_buf, "{s}.pvd", .{self.base_name}) catch
            return error.FilenameTooLong;

        var dir = try std.fs.cwd().openDir(self.output_dir, .{});
        defer dir.close();
        const file = try dir.createFile(pvd_name, .{});
        defer file.close();
        try file.writeAll(output.items);
    }
};

/// Idempotent directory creation: succeeds whether or not the directory
/// already exists. Used by every example before writing snapshots.
pub fn ensureDir(path: []const u8) !void {
    std.fs.cwd().makeDir(path) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const testing = std.testing;

test "Plan disabled when frames is zero" {
    const plan = Plan.fromTotal(100, 0);
    try testing.expectEqual(@as(u32, 0), plan.interval);
    try testing.expectEqual(@as(u32, 0), plan.capacity);
}

test "Plan disabled when total steps is zero" {
    const plan = Plan.fromTotal(0, 5);
    try testing.expectEqual(@as(u32, 0), plan.interval);
    try testing.expectEqual(@as(u32, 0), plan.capacity);
}

test "Plan even division: 100 steps, 10 frames -> interval 10, capacity 11" {
    const plan = Plan.fromTotal(100, 10);
    try testing.expectEqual(@as(u32, 10), plan.interval);
    try testing.expectEqual(@as(u32, 11), plan.capacity);
}

test "Plan rounds interval down for uneven divisions" {
    // 105 / 10 = 10 (integer division). Capacity = 105 / 10 + 1 = 11.
    const plan = Plan.fromTotal(105, 10);
    try testing.expectEqual(@as(u32, 10), plan.interval);
    try testing.expectEqual(@as(u32, 11), plan.capacity);
}

test "Plan with frames > steps yields interval 1" {
    const plan = Plan.fromTotal(5, 10);
    try testing.expectEqual(@as(u32, 1), plan.interval);
    try testing.expectEqual(@as(u32, 6), plan.capacity);
}

test "Plan saturating capacity does not overflow at u32_max" {
    const plan = Plan.fromTotal(std.math.maxInt(u32), 1);
    try testing.expectEqual(@as(u32, std.math.maxInt(u32)), plan.interval);
    try testing.expect(plan.capacity >= 1);
}

test "Series.dueAt returns true exactly at multiples of the interval" {
    var series = try Series.init(testing.allocator, "/tmp", "fixture", Plan.fromTotal(10, 5));
    defer series.deinit();
    // interval = 2, capacity = 6
    try testing.expectEqual(@as(u32, 2), series.plan.interval);
    try testing.expect(!series.dueAt(1));
    try testing.expect(series.dueAt(2));
    try testing.expect(!series.dueAt(3));
    try testing.expect(series.dueAt(4));
    try testing.expect(series.dueAt(10));
}

test "Series.dueAt returns false when disabled" {
    var series = try Series.init(testing.allocator, "/tmp", "fixture", Plan.fromTotal(10, 0));
    defer series.deinit();
    try testing.expect(!series.enabled());
    try testing.expect(!series.dueAt(1));
    try testing.expect(!series.dueAt(10));
}

test "Series.dueAt refuses captures past capacity" {
    // interval = 1, capacity = 4 (3 / 1 + 1)
    var series = try Series.init(testing.allocator, "/tmp", "fixture", Plan.fromTotal(3, 10));
    defer series.deinit();
    series.count = series.plan.capacity;
    try testing.expect(!series.dueAt(1));
}
