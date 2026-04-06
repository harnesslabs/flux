#!/usr/bin/env python3
"""Visualize flux VTK output as a cohesive animation across examples.

Usage:
    uv run tools/visualize.py output/euler_2d
    uv run tools/visualize.py output/euler_2d --field stream_function
    uv run tools/visualize.py output/maxwell_2d --field B_flux --output cavity.png
"""
# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib", "pillow"]
# ///

from __future__ import annotations

import argparse
import glob
import io
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image


plt.rcParams.update(
    {
        "figure.facecolor": "#f7f3ea",
        "axes.facecolor": "#f7f3ea",
        "savefig.facecolor": "#f7f3ea",
        "axes.edgecolor": "#3b352f",
        "axes.labelcolor": "#3b352f",
        "axes.titleweight": "semibold",
        "xtick.color": "#3b352f",
        "ytick.color": "#3b352f",
        "font.size": 10,
    }
)

SCALAR_PREFERENCES = [
    "vorticity",
    "B_flux",
    "E_intensity",
    "stream_function",
]
VECTOR_PREFERENCES = [
    "velocity",
]


@dataclass
class DataArray:
    name: str
    values: np.ndarray
    num_components: int
    association: str


@dataclass
class VtuData:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    dimension: int
    cells: np.ndarray
    point_data: dict[str, DataArray]
    cell_data: dict[str, DataArray]
    timesteps: dict[str, float] | None = None


def parse_vtu(path: str) -> VtuData:
    """Parse a VTU file and return geometry plus point/cell data arrays."""
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError(f"Missing Piece element in {path}")

    n_points = int(piece.attrib["NumberOfPoints"])
    n_cells = int(piece.attrib["NumberOfCells"])

    points_el = piece.find(".//Points/DataArray")
    if points_el is None or points_el.text is None:
        raise ValueError(f"Missing point coordinates in {path}")
    coords = np.fromstring(points_el.text, sep=" ").reshape(n_points, 3)

    connectivity_el = piece.find(".//Cells/DataArray[@Name='connectivity']")
    offsets_el = piece.find(".//Cells/DataArray[@Name='offsets']")
    types_el = piece.find(".//Cells/DataArray[@Name='types']")
    if connectivity_el is None or offsets_el is None or types_el is None:
        raise ValueError(f"Missing cell topology arrays in {path}")

    conn = np.fromstring(connectivity_el.text or "", sep=" ", dtype=int)
    offsets = np.fromstring(offsets_el.text or "", sep=" ", dtype=int)
    cell_types = np.fromstring(types_el.text or "", sep=" ", dtype=int)
    if len(cell_types) == 0:
        raise ValueError(f"No VTK cell types found in {path}")

    vertices_per_cell = offsets[0]
    if np.any(np.diff(offsets) != vertices_per_cell):
        raise ValueError("Mixed cell sizes are not supported")
    cells = conn.reshape(n_cells, vertices_per_cell)

    cell_type = int(cell_types[0])
    if np.any(cell_types != cell_type):
        raise ValueError("Mixed cell types are not supported")

    point_data = parse_data_arrays(piece.findall(".//PointData/DataArray"), n_points, "point")
    cell_data = parse_data_arrays(piece.findall(".//CellData/DataArray"), n_cells, "cell")

    return VtuData(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        dimension=2 if cell_type == 5 else 3,
        cells=cells,
        point_data=point_data,
        cell_data=cell_data,
    )


def parse_data_arrays(elements: list[ET.Element], entity_count: int, association: str) -> dict[str, DataArray]:
    arrays: dict[str, DataArray] = {}
    for da in elements:
        name = da.attrib["Name"]
        num_components = int(da.attrib.get("NumberOfComponents", "1"))
        flat = np.fromstring(da.text or "", sep=" ")
        expected = entity_count * num_components
        if flat.size != expected:
            raise ValueError(
                f"{association} array '{name}' expected {expected} values, found {flat.size}"
            )
        values = flat.reshape(entity_count, num_components) if num_components > 1 else flat
        arrays[name] = DataArray(
            name=name,
            values=values,
            num_components=num_components,
            association=association,
        )
    return arrays


def choose_scalar_field(data: VtuData, requested: str | None) -> DataArray:
    if requested is not None:
        match = find_field(data, requested)
        if match.num_components != 1:
            raise ValueError(f"Field '{requested}' is vector-valued; choose a scalar field")
        return match

    for name in SCALAR_PREFERENCES:
        arrays = [data.cell_data.get(name), data.point_data.get(name)]
        for array in arrays:
            if array is not None and array.num_components == 1:
                return array

    for array in list(data.cell_data.values()) + list(data.point_data.values()):
        if array.num_components == 1:
            return array

    raise ValueError("No scalar field found in VTU output")


def choose_vector_field(data: VtuData, requested: str | None) -> DataArray | None:
    if requested == "none":
        return None
    if requested is not None:
        match = find_field(data, requested)
        if match.num_components <= 1:
            raise ValueError(f"Field '{requested}' is scalar-valued; choose a vector field")
        return match

    for name in VECTOR_PREFERENCES:
        arrays = [data.cell_data.get(name), data.point_data.get(name)]
        for array in arrays:
            if array is not None and array.num_components > 1:
                return array
    return None


def find_field(data: VtuData, name: str) -> DataArray:
    if name in data.cell_data:
        return data.cell_data[name]
    if name in data.point_data:
        return data.point_data[name]
    available = sorted(list(data.cell_data.keys()) + list(data.point_data.keys()))
    raise ValueError(f"Field '{name}' not found. Available: {available}")


def choose_style(field: DataArray) -> tuple[str, bool, str]:
    """Return (cmap, symmetric_range, label)."""
    if field.name in {"vorticity", "B_flux", "stream_function"}:
        return ("RdBu_r", True, field.name.replace("_", " "))
    if field.name in {"E_intensity"}:
        return ("magma", False, field.name.replace("_", " "))
    return ("viridis", False, field.name.replace("_", " "))


def compute_range(values: list[np.ndarray], symmetric: bool) -> tuple[float, float]:
    data_min = min(float(np.min(v)) for v in values)
    data_max = max(float(np.max(v)) for v in values)
    if symmetric:
        bound = max(abs(data_min), abs(data_max))
        if bound == 0.0:
            bound = 1.0
        return -bound, bound
    if data_min == data_max:
        eps = 1.0 if data_min == 0.0 else abs(data_min) * 0.1
        return data_min - eps, data_max + eps
    return data_min, data_max


def cell_barycenters(data: VtuData) -> np.ndarray:
    return np.mean(
        np.stack(
            [
                data.x[data.cells],
                data.y[data.cells],
                data.z[data.cells],
            ],
            axis=-1,
        ),
        axis=1,
    )


def vector_components(field: DataArray) -> tuple[np.ndarray, np.ndarray]:
    if field.association == "cell":
        vectors = field.values[:, :2]
    else:
        vectors = field.values[:, :2]
    return vectors[:, 0], vectors[:, 1]


def render_frame(
    data: VtuData,
    scalar_field: DataArray,
    vector_field: DataArray | None,
    vmin: float,
    vmax: float,
    frame_idx: int,
    n_frames: int,
    series_name: str,
) -> Image.Image:
    cmap, _, colorbar_label = choose_style(scalar_field)
    fig = plt.figure(figsize=(8.2, 6.6), dpi=120)

    if data.dimension == 2:
        ax = fig.add_subplot(1, 1, 1)
        triang = mtri.Triangulation(data.x, data.y, data.cells)
        artist = render_2d_scalar(ax, triang, data, scalar_field, cmap, vmin, vmax)
        overlay_2d_vectors(ax, data, vector_field)
        ax.triplot(triang, color="#ffffff", linewidth=0.18, alpha=0.22)
        ax.set_aspect("equal")
        ax.set_xlim(float(np.min(data.x)), float(np.max(data.x)))
        ax.set_ylim(float(np.min(data.y)), float(np.max(data.y)))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    else:
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        artist = render_3d_scalar(ax, data, scalar_field, cmap, vmin, vmax)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=24, azim=-58)
        ax.set_box_aspect(
            (
                max(np.ptp(data.x), 1e-12),
                max(np.ptp(data.y), 1e-12),
                max(np.ptp(data.z), 1e-12),
            )
        )

    title = f"{series_name}  |  {scalar_field.name.replace('_', ' ')}"
    subtitle = f"frame {frame_idx}/{n_frames}"
    ax.set_title(f"{title}\n{subtitle}", pad=10)

    colorbar = fig.colorbar(artist, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label(colorbar_label)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


def render_2d_scalar(
    ax: plt.Axes,
    triang: mtri.Triangulation,
    data: VtuData,
    field: DataArray,
    cmap: str,
    vmin: float,
    vmax: float,
):
    if field.association == "cell":
        return ax.tripcolor(
            triang,
            field.values,
            shading="flat",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
    return ax.tripcolor(
        triang,
        field.values,
        shading="gouraud",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )


def overlay_2d_vectors(ax: plt.Axes, data: VtuData, field: DataArray | None) -> None:
    if field is None or field.num_components <= 1:
        return

    if field.association == "cell":
        positions = cell_barycenters(data)[:, :2]
    else:
        positions = np.column_stack([data.x, data.y])

    u, v = vector_components(field)
    magnitudes = np.hypot(u, v)
    nonzero = np.flatnonzero(magnitudes > 1e-14)
    if nonzero.size == 0:
        return

    target = min(160, nonzero.size)
    stride = max(1, math.ceil(nonzero.size / target))
    sampled = nonzero[::stride]

    ax.quiver(
        positions[sampled, 0],
        positions[sampled, 1],
        u[sampled],
        v[sampled],
        color="#111111",
        alpha=0.6,
        pivot="mid",
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0022,
    )


def render_3d_scalar(
    ax: plt.Axes,
    data: VtuData,
    field: DataArray,
    cmap: str,
    vmin: float,
    vmax: float,
):
    barycenters = cell_barycenters(data)
    values = field.values if field.association == "cell" else field.values[data.cells].mean(axis=1)
    return ax.scatter(
        barycenters[:, 0],
        barycenters[:, 1],
        barycenters[:, 2],
        c=values,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        s=max(8, int(2400 / max(len(values), 1))),
        depthshade=False,
        edgecolors="none",
    )


def output_format(output_path: str) -> str:
    suffix = Path(output_path).suffix.lower()
    if suffix == ".gif":
        return "gif"
    if suffix in {".png", ".apng"}:
        return "apng"
    raise ValueError(
        f"Unsupported output extension '{suffix}'. Use .gif or .png/.apng."
    )


def save_animation(frames: list[Image.Image], output_path: str, fps: int) -> None:
    duration_ms = int(1000 / fps)
    fmt = output_format(output_path)
    if fmt == "gif":
        frames[0].save(
            output_path,
            save_all=True,
            append_images=frames[1:],
            duration=duration_ms,
            loop=0,
        )
        return

    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        format="PNG",
    )


def detect_series_name(input_dir: str) -> str:
    return Path(input_dir).name.replace("_", " ")


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize flux VTK output as GIF or APNG")
    parser.add_argument("input_dir", help="Directory containing .vtu files")
    parser.add_argument(
        "--field",
        default=None,
        help="Scalar field to plot; default auto-selects a good field for the example",
    )
    parser.add_argument(
        "--vectors",
        default=None,
        help="Optional vector field overlay for 2D plots; default auto-selects velocity when present, use 'none' to disable",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output animation path (default: <input_dir>/animation.png; use .gif for palette-limited GIF)",
    )
    parser.add_argument("--fps", type=int, default=12, help="Frames per second")
    args = parser.parse_args()

    vtu_files = sorted(glob.glob(str(Path(args.input_dir) / "*.vtu")))
    if not vtu_files:
        print(f"No .vtu files found in {args.input_dir}")
        return

    print(f"Found {len(vtu_files)} frames in {args.input_dir}")
    first = parse_vtu(vtu_files[0])
    scalar_field = choose_scalar_field(first, args.field)
    vector_field = choose_vector_field(first, args.vectors)
    _, symmetric_range, _ = choose_style(scalar_field)

    print(f"Scalar field: {scalar_field.name} ({scalar_field.association})")
    if vector_field is not None and first.dimension == 2:
        print(f"Vector overlay: {vector_field.name} ({vector_field.association})")

    scalar_values: list[np.ndarray] = [scalar_field.values]
    parsed_frames: list[VtuData] = [first]
    for path in vtu_files[1:]:
        frame = parse_vtu(path)
        parsed_frames.append(frame)
        scalar_values.append(find_field(frame, scalar_field.name).values)

    vmin, vmax = compute_range(scalar_values, symmetric_range)
    print(f"Value range: [{vmin:.4e}, {vmax:.4e}]")

    series_name = detect_series_name(args.input_dir)
    frames: list[Image.Image] = []
    for i, frame in enumerate(parsed_frames):
        if i % 10 == 0:
            print(f"  rendering frame {i + 1}/{len(parsed_frames)}...")
        scalar = find_field(frame, scalar_field.name)
        vector = None if vector_field is None else find_field(frame, vector_field.name)
        image = render_frame(
            frame,
            scalar,
            vector,
            vmin,
            vmax,
            i + 1,
            len(parsed_frames),
            series_name,
        )
        frames.append(image)

    output_path = args.output or str(Path(args.input_dir) / "animation.png")
    save_animation(frames, output_path, args.fps)
    fmt = output_format(output_path)
    if fmt == "gif":
        print(
            "\nSaved {path} ({frames} frames, {fps} fps, GIF palette up to 256 colors)".format(
                path=output_path,
                frames=len(frames),
                fps=args.fps,
            )
        )
    else:
        print(
            "\nSaved {path} ({frames} frames, {fps} fps, full-color APNG)".format(
                path=output_path,
                frames=len(frames),
                fps=args.fps,
            )
        )


if __name__ == "__main__":
    main()
