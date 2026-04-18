const std = @import("std");

pub const ParseError = error{
    InvalidArgument,
    MissingValue,
};

pub const Parser = struct {
    args: []const [:0]const u8,
    index: usize = 0,

    pub fn init(args: []const [:0]const u8) Parser {
        return .{ .args = args };
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
};

pub inline fn eql(a: []const u8, b: []const u8) bool {
    return std.mem.eql(u8, a, b);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
