//! Shared utilities still used by the fresh example stack.

pub const snapshot = @import("snapshot.zig");
pub const progress = @import("progress.zig");
pub const viz = @import("viz.zig");

pub const Plan = snapshot.Plan;
pub const PlanOptions = snapshot.PlanOptions;
pub const Series = snapshot.Series;
pub const ensureDir = snapshot.ensureDir;
pub const Progress = progress.Progress;
pub const formatDuration = progress.formatDuration;

pub fn framesToInterval(total_steps: u32, frames: u32) u32 {
    if (frames == 0) return total_steps;
    return @max(@as(u32, 1), total_steps / frames);
}

test {
    @import("std").testing.refAllDeclsRecursive(@This());
}
