const vtk = @import("vtk.zig");
const vtk_fields = @import("vtk_fields.zig");
const vtk_series = @import("vtk_series.zig");

pub const DataArraySlice = vtk.DataArraySlice;
pub const parseDataArray = vtk.parseDataArray;
pub const write = vtk.write;

pub const project_edges_to_faces = vtk_fields.project_edges_to_faces;
pub const write_fields = vtk_fields.write_fields;

pub const PvdEntry = vtk_series.PvdEntry;
pub const max_snapshot_filename_length = vtk_series.max_snapshot_filename_length;
pub const snapshot_filename = vtk_series.snapshot_filename;
pub const write_pvd = vtk_series.write_pvd;
