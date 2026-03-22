//! Minimal static file server for Zig autodoc output.
//!
//! Serves `zig-out/docs/` over HTTP so the Wasm-based documentation viewer
//! works in the browser (browsers block Wasm from `file://` due to CORS).
//!
//! Usage: `zig build serve-docs` → opens http://127.0.0.1:8080

const std = @import("std");
const net = std.net;
const mem = std.mem;

const port: u16 = 8080;

const FileEntry = struct {
    path: []const u8,
    content_type: []const u8,
};

const known_files = [_]FileEntry{
    .{ .path = "/", .content_type = "text/html" },
    .{ .path = "/index.html", .content_type = "text/html" },
    .{ .path = "/main.js", .content_type = "application/javascript" },
    .{ .path = "/main.wasm", .content_type = "application/wasm" },
    .{ .path = "/sources.tar", .content_type = "application/x-tar" },
};

pub fn main() !void {
    var gpa_instance: std.heap.GeneralPurposeAllocator(.{}) = .init;
    defer _ = gpa_instance.deinit();
    const gpa = gpa_instance.allocator();

    const docs_dir = std.fs.cwd().openDir("zig-out/docs", .{}) catch |err| {
        std.log.err("Cannot open zig-out/docs/: {s}. Run `zig build docs` first.", .{@errorName(err)});
        std.process.exit(1);
    };
    _ = docs_dir;

    const address = net.Address.parseIp("127.0.0.1", port) catch unreachable;
    var server = try address.listen(.{ .reuse_address = true });

    const url = std.fmt.comptimePrint("http://127.0.0.1:{d}/", .{port});
    std.debug.print("Serving docs at {s}\n", .{url});

    // Open browser
    var child = std.process.Child.init(&.{ "open", url }, gpa);
    _ = child.spawnAndWait() catch {};

    while (true) {
        const connection = server.accept() catch |err| {
            std.log.err("accept failed: {s}", .{@errorName(err)});
            continue;
        };
        _ = std.Thread.spawn(.{}, handleConnection, .{connection}) catch |err| {
            std.log.err("spawn failed: {s}", .{@errorName(err)});
            connection.stream.close();
        };
    }
}

fn handleConnection(connection: net.Server.Connection) void {
    defer connection.stream.close();

    var recv_buffer: [4096]u8 = undefined;
    var send_buffer: [4096]u8 = undefined;
    var conn_reader = connection.stream.reader(&recv_buffer);
    var conn_writer = connection.stream.writer(&send_buffer);
    var server = std.http.Server.init(conn_reader.interface(), &conn_writer.interface);

    while (server.reader.state == .ready) {
        var request = server.receiveHead() catch return;
        serveRequest(&request) catch return;
    }
}

fn serveRequest(request: *std.http.Server.Request) !void {
    const target = request.head.target;

    // Map target to file
    var file_path: []const u8 = undefined;
    var content_type: []const u8 = undefined;
    var found = false;

    for (known_files) |entry| {
        if (mem.eql(u8, target, entry.path)) {
            file_path = if (mem.eql(u8, entry.path, "/")) "index.html" else entry.path[1..];
            content_type = entry.content_type;
            found = true;
            break;
        }
    }

    if (!found) {
        try request.respond("not found", .{
            .status = .not_found,
            .extra_headers = &.{
                .{ .name = "content-type", .value = "text/plain" },
            },
        });
        return;
    }

    // Read file from zig-out/docs/
    var dir = std.fs.cwd().openDir("zig-out/docs", .{}) catch {
        try request.respond("docs not found", .{ .status = .not_found });
        return;
    };
    defer dir.close();

    const contents = dir.readFileAlloc(std.heap.page_allocator, file_path, 50 * 1024 * 1024) catch {
        try request.respond("file not found", .{ .status = .not_found });
        return;
    };
    defer std.heap.page_allocator.free(contents);

    try request.respond(contents, .{
        .extra_headers = &.{
            .{ .name = "content-type", .value = content_type },
            .{ .name = "cache-control", .value = "no-cache" },
        },
    });
}
