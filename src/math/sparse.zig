//! Sparse linear algebra primitives.
//!
//! Provides `CsrMatrix(T)`, a compressed sparse row matrix parameterized on
//! value type. Used for boundary/incidence operators (`CsrMatrix(i8)` with
//! entries in {−1, 0, +1}) and real-valued operators (`CsrMatrix(f64)`).

const std = @import("std");

const testing = std.testing;

/// Compressed sparse row matrix parameterized on value type.
///
/// Standard CSR layout: row pointers index into parallel arrays of column
/// indices and values. Use `CsrMatrix(i8)` for incidence/boundary operators
/// with entries in {−1, 0, +1}, or `CsrMatrix(f64)` for real-valued operators.
pub fn CsrMatrix(comptime T: type) type {
    return struct {
        const Self = @This();

        /// `row_ptr[i]..row_ptr[i+1]` indexes into `col_idx`/`values` for row i.
        /// Length: `n_rows + 1`.
        row_ptr: []u32,
        /// Column indices for each nonzero entry. Length: nnz.
        col_idx: []u32,
        /// Values for each nonzero entry. Length: nnz.
        values: []T,
        n_rows: u32,
        n_cols: u32,

        /// Allocate a CSR matrix with the given dimensions and nonzero capacity.
        /// Caller must populate `row_ptr`, `col_idx`, and `values` after init.
        pub fn init(allocator: std.mem.Allocator, n_rows: u32, n_cols: u32, nonzero_count: u32) !Self {
            return .{
                .row_ptr = try allocator.alloc(u32, @as(usize, n_rows) + 1),
                .col_idx = try allocator.alloc(u32, nonzero_count),
                .values = try allocator.alloc(T, nonzero_count),
                .n_rows = n_rows,
                .n_cols = n_cols,
            };
        }

        /// Free all allocated storage.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.row_ptr);
            allocator.free(self.col_idx);
            allocator.free(self.values);
        }

        /// Number of nonzero entries.
        ///
        /// Asserts that the count fits in u32. This is always true for matrices
        /// constructed via `init` (which takes a u32 count), but protects against
        /// a corrupted or externally-constructed matrix.
        pub fn nnz(self: Self) u32 {
            std.debug.assert(self.col_idx.len <= std.math.maxInt(u32));
            return @intCast(self.col_idx.len);
        }

        /// Column indices and values for a given row.
        pub fn row(self: Self, r: u32) struct { cols: []const u32, vals: []const T } {
            const start = self.row_ptr[r];
            const end = self.row_ptr[r + 1];
            return .{
                .cols = self.col_idx[start..end],
                .vals = self.values[start..end],
            };
        }

        /// Compute y += Aᵀ x (transpose sparse matrix–vector product).
        ///
        /// Iterates over rows of A, scattering each row's contribution into
        /// the output vector. Caller must zero-initialize `output` before
        /// calling. `input` has length `n_rows`, `output` has length `n_cols`.
        pub fn transpose_multiply(self: Self, input_vals: []const f64, output: []f64) void {
            std.debug.assert(input_vals.len == self.n_rows);
            std.debug.assert(output.len == self.n_cols);
            for (0..self.n_rows) |row_idx| {
                const r = self.row(@intCast(row_idx));
                for (r.cols, r.vals) |col, val| {
                    output[col] += @as(f64, @floatFromInt(val)) * input_vals[row_idx];
                }
            }
        }
    };
}

/// Compute y = A x (forward sparse matrix–vector product) for f64 matrices.
///
/// Standard CSR SpMV: for each row, dot product of the row with the
/// input vector. `input` has length `n_cols`, `output` has length `n_rows`.
/// Output is overwritten (not accumulated).
pub fn spmv(matrix: CsrMatrix(f64), input_vals: []const f64, output: []f64) void {
    std.debug.assert(input_vals.len == matrix.n_cols);
    std.debug.assert(output.len == matrix.n_rows);
    for (0..matrix.n_rows) |row_idx| {
        const r = matrix.row(@intCast(row_idx));
        var sum: f64 = 0;
        for (r.cols, r.vals) |col, val| {
            sum += val * input_vals[col];
        }
        output[row_idx] = sum;
    }
}

/// Packed sign stream for incidence-valued sparse matrices.
///
/// CSR already omits zeros, so incidence values only need one bit per stored
/// entry: 0 => -1, 1 => +1.
pub const PackedIncidenceSigns = struct {
    bytes: []u8,
    len: u32,

    pub fn bytesForEntries(entry_count: u32) usize {
        return @intCast(@divFloor(@as(u64, entry_count) + 7, 8));
    }

    pub fn init(allocator: std.mem.Allocator, entry_count: u32) !PackedIncidenceSigns {
        const byte_count = bytesForEntries(entry_count);
        const bytes = try allocator.alloc(u8, byte_count);
        @memset(bytes, 0);
        return .{
            .bytes = bytes,
            .len = entry_count,
        };
    }

    pub fn deinit(self: *PackedIncidenceSigns, allocator: std.mem.Allocator) void {
        allocator.free(self.bytes);
    }

    pub fn fromSigns(allocator: std.mem.Allocator, signs: []const i8) !PackedIncidenceSigns {
        std.debug.assert(signs.len <= std.math.maxInt(u32));

        var packed_signs = try init(allocator, @intCast(signs.len));
        errdefer packed_signs.deinit(allocator);

        for (signs, 0..) |sign, entry_idx| {
            packed_signs.set(@intCast(entry_idx), sign);
        }

        return packed_signs;
    }

    pub fn get(self: PackedIncidenceSigns, entry_idx: u32) i8 {
        std.debug.assert(entry_idx < self.len);
        const mask = bitMask(entry_idx);
        const byte = self.bytes[byteIndex(entry_idx)];
        return if ((byte & mask) == 0) -1 else 1;
    }

    fn set(self: PackedIncidenceSigns, entry_idx: u32, sign: i8) void {
        std.debug.assert(entry_idx < self.len);
        std.debug.assert(sign == -1 or sign == 1);

        const byte_idx = byteIndex(entry_idx);
        const mask = bitMask(entry_idx);
        if (sign == 1) {
            self.bytes[byte_idx] |= mask;
        } else {
            self.bytes[byte_idx] &= ~mask;
        }
    }

    fn byteIndex(entry_idx: u32) usize {
        return @intCast(@divFloor(entry_idx, 8));
    }

    fn bitMask(entry_idx: u32) u8 {
        const bit_idx: u3 = @intCast(@mod(entry_idx, 8));
        return @as(u8, 1) << bit_idx;
    }
};

/// Prototype CSR variant for incidence matrices with packed sign storage.
pub const PackedIncidenceMatrix = struct {
    pub const Row = struct {
        cols: []const u32,
        start: u32,
        signs: PackedIncidenceSigns,

        pub fn sign(self: Row, local_idx: usize) i8 {
            std.debug.assert(local_idx < self.cols.len);
            return self.signs.get(self.start + @as(u32, @intCast(local_idx)));
        }
    };

    row_ptr: []u32,
    col_idx: []u32,
    signs: PackedIncidenceSigns,
    n_rows: u32,
    n_cols: u32,

    pub fn fromCsr(allocator: std.mem.Allocator, matrix: CsrMatrix(i8)) !PackedIncidenceMatrix {
        const row_ptr = try allocator.dupe(u32, matrix.row_ptr);
        errdefer allocator.free(row_ptr);

        const col_idx = try allocator.dupe(u32, matrix.col_idx);
        errdefer allocator.free(col_idx);

        for (matrix.values) |value| {
            std.debug.assert(value == -1 or value == 1);
        }

        const signs = try PackedIncidenceSigns.fromSigns(allocator, matrix.values);
        errdefer {
            var signs_copy = signs;
            signs_copy.deinit(allocator);
        }

        return .{
            .row_ptr = row_ptr,
            .col_idx = col_idx,
            .signs = signs,
            .n_rows = matrix.n_rows,
            .n_cols = matrix.n_cols,
        };
    }

    pub fn deinit(self: *PackedIncidenceMatrix, allocator: std.mem.Allocator) void {
        allocator.free(self.row_ptr);
        allocator.free(self.col_idx);
        self.signs.deinit(allocator);
    }

    pub fn row(self: PackedIncidenceMatrix, r: u32) Row {
        const start = self.row_ptr[r];
        const end = self.row_ptr[r + 1];
        return .{
            .cols = self.col_idx[start..end],
            .start = start,
            .signs = self.signs,
        };
    }

    pub fn transpose_multiply(self: PackedIncidenceMatrix, input_vals: []const f64, output: []f64) void {
        std.debug.assert(input_vals.len == self.n_rows);
        std.debug.assert(output.len == self.n_cols);

        for (0..self.n_rows) |row_idx_usize| {
            const row_idx: usize = row_idx_usize;
            const row_view = self.row(@intCast(row_idx_usize));
            const input_value = input_vals[row_idx];

            for (row_view.cols, 0..) |col, entry_offset| {
                const sign = row_view.sign(entry_offset);
                output[col] += if (sign == 1) input_value else -input_value;
            }
        }
    }
};

/// COO (coordinate) format assembler for building CSR matrices incrementally.
///
/// Use `addEntry` to accumulate element contributions in any order. Duplicate
/// (row, col) entries are summed during `build`, which produces a sorted CSR
/// matrix. Designed for finite element assembly where each element contributes
/// to a small block of the global matrix.
pub fn TripletAssembler(comptime T: type) type {
    return struct {
        const Self = @This();

        const Triplet = struct {
            row: u32,
            col: u32,
            val: T,
        };

        triplets: std.ArrayListUnmanaged(Triplet),
        n_rows: u32,
        n_cols: u32,

        pub fn init(n_rows: u32, n_cols: u32) Self {
            return .{
                .triplets = .{},
                .n_rows = n_rows,
                .n_cols = n_cols,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.triplets.deinit(allocator);
        }

        /// Add a single (row, col, value) entry. Duplicates are summed in build().
        pub fn addEntry(self: *Self, allocator: std.mem.Allocator, row_idx: u32, col: u32, val: T) !void {
            std.debug.assert(row_idx < self.n_rows);
            std.debug.assert(col < self.n_cols);
            try self.triplets.append(allocator, .{ .row = row_idx, .col = col, .val = val });
        }

        /// Build a CSR matrix from the accumulated triplets.
        ///
        /// Sorts by (row, col), sums duplicates, and produces a compact CSR.
        /// The assembler can still be used after build (e.g., to build a
        /// second matrix with different values but the same pattern).
        pub fn build(self: *Self, allocator: std.mem.Allocator) !CsrMatrix(T) {
            const items = self.triplets.items;

            // Sort by (row, col) for CSR construction.
            std.mem.sort(Triplet, items, {}, struct {
                fn lessThan(_: void, a: Triplet, b: Triplet) bool {
                    if (a.row != b.row) return a.row < b.row;
                    return a.col < b.col;
                }
            }.lessThan);

            // Count unique (row, col) pairs to determine nnz.
            var unique_count: u32 = 0;
            for (items, 0..) |entry, i| {
                if (i == 0 or entry.row != items[i - 1].row or entry.col != items[i - 1].col) {
                    unique_count += 1;
                }
            }

            var matrix = try CsrMatrix(T).init(allocator, self.n_rows, self.n_cols, unique_count);
            errdefer matrix.deinit(allocator);

            // Fill CSR arrays — sum duplicates.
            var write_idx: u32 = 0;
            var current_row: u32 = 0;
            matrix.row_ptr[0] = 0;

            for (items, 0..) |entry, i| {
                // Advance row pointers for any skipped rows.
                while (current_row < entry.row) {
                    current_row += 1;
                    matrix.row_ptr[current_row] = write_idx;
                }

                if (i > 0 and entry.row == items[i - 1].row and entry.col == items[i - 1].col) {
                    // Duplicate — sum into existing entry.
                    matrix.values[write_idx - 1] += entry.val;
                } else {
                    // New unique entry.
                    matrix.col_idx[write_idx] = entry.col;
                    matrix.values[write_idx] = entry.val;
                    write_idx += 1;
                }
            }

            // Fill remaining row pointers.
            while (current_row < self.n_rows) {
                current_row += 1;
                matrix.row_ptr[current_row] = write_idx;
            }

            std.debug.assert(write_idx == unique_count);
            return matrix;
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

test "CsrMatrix init and nnz" {
    const allocator = testing.allocator;
    var m = try CsrMatrix(f64).init(allocator, 3, 4, 5);
    defer m.deinit(allocator);

    try testing.expectEqual(@as(u32, 3), m.n_rows);
    try testing.expectEqual(@as(u32, 4), m.n_cols);
    try testing.expectEqual(@as(u32, 5), m.nnz());
}

test "CsrMatrix row returns correct slice" {
    const allocator = testing.allocator;
    // 2×3 matrix: row 0 has 2 entries, row 1 has 1 entry
    var m = try CsrMatrix(i8).init(allocator, 2, 3, 3);
    defer m.deinit(allocator);

    m.row_ptr[0] = 0;
    m.row_ptr[1] = 2;
    m.row_ptr[2] = 3;
    m.col_idx[0] = 0;
    m.col_idx[1] = 2;
    m.col_idx[2] = 1;
    m.values[0] = 1;
    m.values[1] = -1;
    m.values[2] = 1;

    const r0 = m.row(0);
    try testing.expectEqual(@as(usize, 2), r0.cols.len);
    try testing.expectEqual(@as(u32, 0), r0.cols[0]);
    try testing.expectEqual(@as(u32, 2), r0.cols[1]);
    try testing.expectEqual(@as(i8, 1), r0.vals[0]);
    try testing.expectEqual(@as(i8, -1), r0.vals[1]);

    const r1 = m.row(1);
    try testing.expectEqual(@as(usize, 1), r1.cols.len);
    try testing.expectEqual(@as(u32, 1), r1.cols[0]);
    try testing.expectEqual(@as(i8, 1), r1.vals[0]);
}

test "spmv (forward SpMV)" {
    const allocator = testing.allocator;
    // 2×3 matrix: [[1.0, 0.0, -1.0], [0.0, 2.0, 0.0]]
    var m = try CsrMatrix(f64).init(allocator, 2, 3, 3);
    defer m.deinit(allocator);

    m.row_ptr[0] = 0;
    m.row_ptr[1] = 2;
    m.row_ptr[2] = 3;
    m.col_idx[0] = 0;
    m.col_idx[1] = 2;
    m.col_idx[2] = 1;
    m.values[0] = 1.0;
    m.values[1] = -1.0;
    m.values[2] = 2.0;

    // y = A x, where x = [3.0, 5.0, -2.0]
    // y[0] = 1.0*3.0 + (-1.0)*(-2.0) = 5.0
    // y[1] = 2.0*5.0 = 10.0
    const input = [_]f64{ 3.0, 5.0, -2.0 };
    var output = [_]f64{ 0.0, 0.0 };
    spmv(m, &input, &output);

    try testing.expectApproxEqAbs(@as(f64, 5.0), output[0], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 10.0), output[1], 1e-15);
}

test "CsrMatrix transpose_multiply" {
    const allocator = testing.allocator;
    // 2×3 matrix: [[1, 0, -1], [0, 2, 0]]
    var m = try CsrMatrix(i8).init(allocator, 2, 3, 3);
    defer m.deinit(allocator);

    m.row_ptr[0] = 0;
    m.row_ptr[1] = 2;
    m.row_ptr[2] = 3;
    m.col_idx[0] = 0;
    m.col_idx[1] = 2;
    m.col_idx[2] = 1;
    m.values[0] = 1;
    m.values[1] = -1;
    m.values[2] = 2;

    // y = Aᵀ x, where x = [3.0, 5.0]
    // Aᵀ = [[1, 0], [0, 2], [-1, 0]]
    // y = [3.0, 10.0, -3.0]
    const input = [_]f64{ 3.0, 5.0 };
    var output = [_]f64{ 0.0, 0.0, 0.0 };
    m.transpose_multiply(&input, &output);

    try testing.expectApproxEqAbs(@as(f64, 3.0), output[0], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 10.0), output[1], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, -3.0), output[2], 1e-15);
}

// ═══════════════════════════════════════════════════════════════════════════
// TripletAssembler tests
// ═══════════════════════════════════════════════════════════════════════════

test "TripletAssembler builds correct CSR from unsorted entries" {
    const allocator = testing.allocator;
    // Build a 3×3 matrix: [[1, 0, 2], [0, 3, 0], [4, 0, 5]]
    var assembler = TripletAssembler(f64).init(3, 3);
    defer assembler.deinit(allocator);

    // Add entries out of order.
    try assembler.addEntry(allocator, 2, 2, 5.0);
    try assembler.addEntry(allocator, 0, 0, 1.0);
    try assembler.addEntry(allocator, 1, 1, 3.0);
    try assembler.addEntry(allocator, 2, 0, 4.0);
    try assembler.addEntry(allocator, 0, 2, 2.0);

    var m = try assembler.build(allocator);
    defer m.deinit(allocator);

    try testing.expectEqual(@as(u32, 3), m.n_rows);
    try testing.expectEqual(@as(u32, 3), m.n_cols);
    try testing.expectEqual(@as(u32, 5), m.nnz());

    // Verify via SpMV: y = A * [1, 1, 1] = [3, 3, 9]
    const x = [_]f64{ 1.0, 1.0, 1.0 };
    var y = [_]f64{ 0.0, 0.0, 0.0 };
    spmv(m, &x, &y);

    try testing.expectApproxEqAbs(@as(f64, 3.0), y[0], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 3.0), y[1], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 9.0), y[2], 1e-15);
}

test "TripletAssembler sums duplicate entries" {
    const allocator = testing.allocator;
    var assembler = TripletAssembler(f64).init(2, 2);
    defer assembler.deinit(allocator);

    // Two contributions to (0,0) and (1,1) — simulates element assembly.
    try assembler.addEntry(allocator, 0, 0, 1.0);
    try assembler.addEntry(allocator, 0, 0, 2.0);
    try assembler.addEntry(allocator, 1, 1, 3.0);
    try assembler.addEntry(allocator, 1, 1, 4.0);
    try assembler.addEntry(allocator, 0, 1, 0.5);

    var m = try assembler.build(allocator);
    defer m.deinit(allocator);

    // Only 3 unique entries after summing duplicates.
    try testing.expectEqual(@as(u32, 3), m.nnz());

    // (0,0) = 1+2 = 3, (0,1) = 0.5, (1,1) = 3+4 = 7
    const x = [_]f64{ 1.0, 1.0 };
    var y = [_]f64{ 0.0, 0.0 };
    spmv(m, &x, &y);

    try testing.expectApproxEqAbs(@as(f64, 3.5), y[0], 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 7.0), y[1], 1e-15);
}

test "TripletAssembler handles empty rows" {
    const allocator = testing.allocator;
    var assembler = TripletAssembler(f64).init(3, 3);
    defer assembler.deinit(allocator);

    // Only row 0 and row 2 have entries — row 1 is empty.
    try assembler.addEntry(allocator, 0, 0, 1.0);
    try assembler.addEntry(allocator, 2, 2, 2.0);

    var m = try assembler.build(allocator);
    defer m.deinit(allocator);

    try testing.expectEqual(@as(u32, 2), m.nnz());

    // Row 1 should have no entries.
    const r1 = m.row(1);
    try testing.expectEqual(@as(usize, 0), r1.cols.len);
}

test "PackedIncidenceSigns round-trips random signs" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x51A9_B17);

    const lengths = [_]u32{ 1, 2, 7, 8, 9, 31, 32, 33, 257 };
    for (lengths) |len| {
        const signs = try allocator.alloc(i8, len);
        defer allocator.free(signs);

        for (signs) |*sign| {
            sign.* = if (rng.random().boolean()) @as(i8, 1) else @as(i8, -1);
        }

        var packed_signs = try PackedIncidenceSigns.fromSigns(allocator, signs);
        defer packed_signs.deinit(allocator);

        try testing.expectEqual(PackedIncidenceSigns.bytesForEntries(len), packed_signs.bytes.len);
        for (signs, 0..) |expected, entry_idx| {
            try testing.expectEqual(expected, packed_signs.get(@intCast(entry_idx)));
        }
    }
}

test "PackedIncidenceMatrix transpose_multiply matches i8 CSR baseline" {
    const allocator = testing.allocator;
    var rng = std.Random.DefaultPrng.init(0x1A61_D3CE);

    inline for (0..32) |_| {
        const n_rows: u32 = 12;
        const n_cols: u32 = 19;
        const max_row_nnz: u32 = 4;

        var matrix = try CsrMatrix(i8).init(allocator, n_rows, n_cols, n_rows * max_row_nnz);
        defer matrix.deinit(allocator);

        var write_idx: u32 = 0;
        matrix.row_ptr[0] = 0;
        for (0..n_rows) |row_idx_usize| {
            const row_idx: u32 = @intCast(row_idx_usize);
            const row_nnz = rng.random().uintLessThanBiased(u32, max_row_nnz + 1);
            var used_cols = [_]bool{false} ** n_cols;
            var selected: u32 = 0;
            while (selected < row_nnz) {
                const col = rng.random().uintLessThanBiased(u32, n_cols);
                if (used_cols[col]) continue;
                used_cols[col] = true;
                selected += 1;
            }

            for (used_cols, 0..) |present, col| {
                if (!present) continue;
                matrix.col_idx[write_idx] = @intCast(col);
                matrix.values[write_idx] = if (rng.random().boolean()) @as(i8, 1) else @as(i8, -1);
                write_idx += 1;
            }
            matrix.row_ptr[row_idx + 1] = write_idx;
        }

        const matrix_view = CsrMatrix(i8){
            .row_ptr = matrix.row_ptr,
            .col_idx = matrix.col_idx[0..write_idx],
            .values = matrix.values[0..write_idx],
            .n_rows = matrix.n_rows,
            .n_cols = matrix.n_cols,
        };

        var packed_matrix = try PackedIncidenceMatrix.fromCsr(allocator, matrix_view);
        defer packed_matrix.deinit(allocator);

        const input = try allocator.alloc(f64, n_rows);
        defer allocator.free(input);
        const dense_output = try allocator.alloc(f64, n_cols);
        defer allocator.free(dense_output);
        const packed_output = try allocator.alloc(f64, n_cols);
        defer allocator.free(packed_output);

        for (input) |*value| {
            value.* = rng.random().float(f64) * 20.0 - 10.0;
        }

        @memset(dense_output, 0.0);
        @memset(packed_output, 0.0);
        matrix_view.transpose_multiply(input, dense_output);
        packed_matrix.transpose_multiply(input, packed_output);

        for (dense_output, packed_output) |expected, actual| {
            try testing.expectApproxEqAbs(expected, actual, 1e-12);
        }
    }
}
