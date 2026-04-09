//! Shared command-line parser for the flux example suite.
//!
//! Every example accepts the same core flags (`--steps`, `--dt`, `--output`,
//! `--frames`, `--grid`, `--domain`, `--refinement`, `--final-time`) so users
//! do not have to relearn the schema per simulation. Each subcommand layers
//! its own physics-specific flags on top by calling `tryCommon` first and
//! falling through on unrecognized arguments.
//!
//! All numeric flags are mandatory-value: `--steps 100`, never `--steps=100`.
//! The two forms could be supported but the explicit space variant matches
//! the existing convention across the project and avoids parser ambiguity.

const std = @import("std");

pub const ParseError = error{
    InvalidArgument,
    MissingValue,
    UnknownFlag,
};

/// Subset of CLI options shared by every example. Fields are optional so the
/// example can distinguish "user did not pass this flag" from "user passed
/// the default value", which matters because some flags have physics-derived
/// defaults that should only be overridden when explicitly requested (notably
/// `--dt`, where the timestep is normally computed from the Courant condition).
pub const Common = struct {
    dt: ?f64 = null,
    steps: ?u32 = null,
    output_dir: ?[]const u8 = null,
    frames: ?u32 = null,
    grid: ?u32 = null,
    domain: ?f64 = null,
    refinement: ?u32 = null,
    final_time: ?f64 = null,
    help: bool = false,
};

pub fn applySharedFields(cfg: anytype, common: Common) void {
    if (common.steps) |value| assignIfPresent(cfg, "steps", value);
    if (common.output_dir) |value| assignIfPresent(cfg, "output_dir", value);
    if (common.frames) |value| assignIfPresent(cfg, "frames", value);
    if (common.grid) |value| assignIfPresent(cfg, "grid", value);
    if (common.domain) |value| assignIfPresent(cfg, "domain", value);
    if (common.refinement) |value| assignIfPresent(cfg, "refinement", value);
    if (common.final_time) |value| assignIfPresent(cfg, "final_time", value);
}

pub fn framesToInterval(total_steps: u32, frames: u32) u32 {
    if (frames == 0) return 0;
    return @max(@as(u32, 1), total_steps / frames);
}

fn assignIfPresent(cfg: anytype, comptime field_name: []const u8, value: anytype) void {
    const ConfigType = @TypeOf(cfg.*);
    if (!@hasField(ConfigType, field_name)) return;

    const FieldType = @FieldType(ConfigType, field_name);
    const ValueType = @TypeOf(value);
    if (FieldType == ValueType) {
        @field(cfg, field_name) = value;
        return;
    }
    if (FieldType == ?ValueType) {
        @field(cfg, field_name) = value;
        return;
    }
}

/// Stateful argument cursor. Holds a slice of argv-style strings and an
/// index. The parser does not own the storage — the caller (typically
/// `std.process.argsAlloc`) is responsible for lifetime.
pub const Parser = struct {
    args: []const [:0]const u8,
    index: usize,

    pub fn init(args: []const [:0]const u8) Parser {
        return .{ .args = args, .index = 0 };
    }

    pub fn remaining(self: Parser) usize {
        return self.args.len - self.index;
    }

    pub fn peek(self: Parser) ?[]const u8 {
        if (self.index >= self.args.len) return null;
        return self.args[self.index];
    }

    pub fn next(self: *Parser) ?[]const u8 {
        if (self.index >= self.args.len) return null;
        const value = self.args[self.index];
        self.index += 1;
        return value;
    }

    pub fn requireValue(self: *Parser, flag: []const u8) ParseError![]const u8 {
        return self.next() orelse {
            std.debug.print("error: {s} requires a value\n", .{flag});
            return ParseError.MissingValue;
        };
    }

    pub fn parseU32(self: *Parser, flag: []const u8) ParseError!u32 {
        const value = try self.requireValue(flag);
        return std.fmt.parseInt(u32, value, 10) catch {
            std.debug.print("error: invalid {s} value '{s}'\n", .{ flag, value });
            return ParseError.InvalidArgument;
        };
    }

    pub fn parseF64(self: *Parser, flag: []const u8) ParseError!f64 {
        const value = try self.requireValue(flag);
        return std.fmt.parseFloat(f64, value) catch {
            std.debug.print("error: invalid {s} value '{s}'\n", .{ flag, value });
            return ParseError.InvalidArgument;
        };
    }

    /// If `arg` matches a known common flag, consume its value (if any) and
    /// store it in `common`. Returns true on a successful match, false when
    /// the flag is unrecognized so the caller can dispatch to its own
    /// example-specific handler. `--help` and `-h` are detected eagerly here
    /// but only set the `help` boolean — the example controls how usage is
    /// printed.
    pub fn tryCommon(self: *Parser, arg: []const u8, common: *Common) ParseError!bool {
        if (eql(arg, "--help") or eql(arg, "-h")) {
            common.help = true;
            return true;
        }
        if (eql(arg, "--dt")) {
            common.dt = try self.parseF64("--dt");
            return true;
        }
        if (eql(arg, "--steps")) {
            common.steps = try self.parseU32("--steps");
            return true;
        }
        if (eql(arg, "--output")) {
            common.output_dir = try self.requireValue("--output");
            return true;
        }
        if (eql(arg, "--frames")) {
            common.frames = try self.parseU32("--frames");
            return true;
        }
        if (eql(arg, "--grid")) {
            common.grid = try self.parseU32("--grid");
            return true;
        }
        if (eql(arg, "--domain")) {
            common.domain = try self.parseF64("--domain");
            return true;
        }
        if (eql(arg, "--refinement")) {
            common.refinement = try self.parseU32("--refinement");
            return true;
        }
        if (eql(arg, "--final-time")) {
            common.final_time = try self.parseF64("--final-time");
            return true;
        }
        return false;
    }
};

inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

const testing = std.testing;

fn makeArgs(comptime literals: []const [:0]const u8) []const [:0]const u8 {
    return literals;
}

test "parser yields nothing on empty args" {
    var parser = Parser.init(&.{});
    try testing.expectEqual(@as(?[]const u8, null), parser.peek());
    try testing.expectEqual(@as(?[]const u8, null), parser.next());
    try testing.expectEqual(@as(usize, 0), parser.remaining());
}

test "parser advances index on next" {
    const args = makeArgs(&.{ "alpha", "beta", "gamma" });
    var parser = Parser.init(args);
    try testing.expectEqualStrings("alpha", parser.next().?);
    try testing.expectEqualStrings("beta", parser.next().?);
    try testing.expectEqualStrings("gamma", parser.next().?);
    try testing.expectEqual(@as(?[]const u8, null), parser.next());
}

test "parseU32 reads following argument" {
    const args = makeArgs(&.{"42"});
    var parser = Parser.init(args);
    try testing.expectEqual(@as(u32, 42), try parser.parseU32("--n"));
}

test "parseU32 rejects non-numeric value" {
    const args = makeArgs(&.{"abc"});
    var parser = Parser.init(args);
    try testing.expectError(ParseError.InvalidArgument, parser.parseU32("--n"));
}

test "parseU32 reports missing value" {
    const args = makeArgs(&.{});
    var parser = Parser.init(args);
    try testing.expectError(ParseError.MissingValue, parser.parseU32("--n"));
}

test "parseF64 reads scientific notation" {
    const args = makeArgs(&.{"1.5e-3"});
    var parser = Parser.init(args);
    try testing.expectEqual(@as(f64, 1.5e-3), try parser.parseF64("--dt"));
}

test "tryCommon recognizes every shared flag" {
    const args = makeArgs(&.{
        "--dt",         "0.01",
        "--steps",      "100",
        "--output",     "tmp/out",
        "--frames",     "10",
        "--grid",       "32",
        "--domain",     "2.0",
        "--refinement", "3",
        "--final-time", "0.5",
    });
    var parser = Parser.init(args);
    var common = Common{};
    while (parser.next()) |arg| {
        try testing.expect(try parser.tryCommon(arg, &common));
    }
    try testing.expectEqual(@as(?f64, 0.01), common.dt);
    try testing.expectEqual(@as(?u32, 100), common.steps);
    try testing.expectEqualStrings("tmp/out", common.output_dir.?);
    try testing.expectEqual(@as(?u32, 10), common.frames);
    try testing.expectEqual(@as(?u32, 32), common.grid);
    try testing.expectEqual(@as(?f64, 2.0), common.domain);
    try testing.expectEqual(@as(?u32, 3), common.refinement);
    try testing.expectEqual(@as(?f64, 0.5), common.final_time);
    try testing.expect(!common.help);
}

test "tryCommon detects --help and -h" {
    {
        const args = makeArgs(&.{"--help"});
        var parser = Parser.init(args);
        var common = Common{};
        try testing.expect(try parser.tryCommon(parser.next().?, &common));
        try testing.expect(common.help);
    }
    {
        const args = makeArgs(&.{"-h"});
        var parser = Parser.init(args);
        var common = Common{};
        try testing.expect(try parser.tryCommon(parser.next().?, &common));
        try testing.expect(common.help);
    }
}

test "tryCommon returns false on unknown flag" {
    const args = makeArgs(&.{"--unknown"});
    var parser = Parser.init(args);
    var common = Common{};
    try testing.expect(!try parser.tryCommon(parser.next().?, &common));
}

test "tryCommon leaves untouched fields null" {
    const args = makeArgs(&.{ "--steps", "5" });
    var parser = Parser.init(args);
    var common = Common{};
    _ = try parser.tryCommon(parser.next().?, &common);
    try testing.expectEqual(@as(?u32, 5), common.steps);
    try testing.expectEqual(@as(?f64, null), common.dt);
    try testing.expectEqual(@as(?u32, null), common.grid);
}
