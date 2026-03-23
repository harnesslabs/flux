//! Sparse linear algebra primitives.
//!
//! Provides `CsrMatrix(T)`, a compressed sparse row matrix parameterized on
//! value type. Used for boundary/incidence operators (`CsrMatrix(i8)` with
//! entries in {−1, 0, +1}) and real-valued operators (`CsrMatrix(f64)`).

const std = @import("std");

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
        pub fn nnz(self: Self) u32 {
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
