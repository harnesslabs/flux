#!/usr/bin/env python3
"""Visualize flux VTK output as an animation.

Usage:
    uv run tools/visualize.py /tmp/flux-cavity
    uv run tools/visualize.py /tmp/flux-dipole --field B_flux
    uv run tools/visualize.py /tmp/flux-cavity --output cavity.gif --fps 15
    uv run tools/visualize.py /tmp/flux-cavity --output cavity.png
"""
# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib", "pillow"]
# ///

import argparse
import glob
import io
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def parse_vtu(path: str) -> dict:
    """Parse a VTK .vtu XML file and return points, cells, and cell data."""
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")

    n_points = int(piece.attrib["NumberOfPoints"])
    n_cells = int(piece.attrib["NumberOfCells"])

    # Points (vertex coordinates).
    points_el = piece.find(".//Points/DataArray")
    coords = np.fromstring(points_el.text, sep=" ").reshape(n_points, 3)

    # Cell connectivity.
    connectivity_el = piece.find(".//Cells/DataArray[@Name='connectivity']")
    conn = np.fromstring(connectivity_el.text, sep=" ", dtype=int)

    offsets_el = piece.find(".//Cells/DataArray[@Name='offsets']")
    offsets = np.fromstring(offsets_el.text, sep=" ", dtype=int)

    types_el = piece.find(".//Cells/DataArray[@Name='types']")
    cell_types = np.fromstring(types_el.text, sep=" ", dtype=int)
    if len(cell_types) == 0:
        raise ValueError(f"No VTK cell types found in {path}")

    vertices_per_cell = offsets[0]
    if np.any(np.diff(offsets) != vertices_per_cell):
        raise ValueError("Mixed cell sizes are not supported")

    cells = conn.reshape(n_cells, vertices_per_cell)
    cell_type = int(cell_types[0])
    if np.any(cell_types != cell_type):
        raise ValueError("Mixed cell types are not supported")

    # Cell data arrays.
    cell_data = {}
    for da in piece.findall(".//CellData/DataArray"):
        name = da.attrib["Name"]
        cell_data[name] = np.fromstring(da.text, sep=" ")

    return {
        "x": coords[:, 0],
        "y": coords[:, 1],
        "z": coords[:, 2],
        "dimension": 2 if cell_type == 5 else 3,
        "cells": cells,
        "cell_data": cell_data,
    }


def render_frame(data: dict, field_name: str, vmin: float, vmax: float,
                 frame_idx: int, n_frames: int) -> Image.Image:
    """Render a single frame as a PIL Image."""
    values = data["cell_data"][field_name]
    fig = plt.figure(figsize=(6, 5), dpi=100)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))

    if data["dimension"] == 2:
        ax = fig.add_subplot(1, 1, 1)
        triang = tri.Triangulation(data["x"], data["y"], data["cells"])
        ax.tripcolor(triang, values, shading="flat", cmap="RdBu_r",
                     vmin=vmin, vmax=vmax)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    else:
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        barycenters = np.mean(
            np.stack(
                [
                    data["x"][data["cells"]],
                    data["y"][data["cells"]],
                    data["z"][data["cells"]],
                ],
                axis=-1,
            ),
            axis=1,
        )
        scatter = ax.scatter(
            barycenters[:, 0],
            barycenters[:, 1],
            barycenters[:, 2],
            c=values,
            cmap="RdBu_r",
            vmin=vmin,
            vmax=vmax,
            s=max(8, int(2400 / max(len(values), 1))),
            depthshade=False,
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=24, azim=-58)
        ax.set_box_aspect((
            max(np.ptp(data["x"]), 1e-12),
            max(np.ptp(data["y"]), 1e-12),
            max(np.ptp(data["z"]), 1e-12),
        ))
        sm = scatter

    ax.set_title(f"{field_name}  (frame {frame_idx}/{n_frames})")
    fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)

    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


def output_format(output_path: str) -> str:
    """Infer animation format from the output filename."""
    suffix = Path(output_path).suffix.lower()
    if suffix == ".gif":
        return "gif"
    if suffix in {".png", ".apng"}:
        return "apng"
    raise ValueError(
        f"Unsupported output extension '{suffix}'. Use .gif for GIF or .png/.apng for APNG."
    )


def save_animation(frames: list[Image.Image], output_path: str, fps: int) -> None:
    """Save frames as either GIF or APNG."""
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


def main():
    parser = argparse.ArgumentParser(description="Visualize flux VTK output as GIF or APNG")
    parser.add_argument("input_dir", help="Directory containing .vtu files")
    parser.add_argument("--field", default="B_flux",
                        help="CellData field to plot (default: B_flux)")
    parser.add_argument("--output", default=None,
                        help="Output animation path (default: <input_dir>/animation.png; use .gif for palette-limited GIF)")
    parser.add_argument("--fps", type=int, default=12, help="Frames per second")
    args = parser.parse_args()

    # Find all .vtu files sorted by name.
    vtu_files = sorted(glob.glob(str(Path(args.input_dir) / "*.vtu")))
    if not vtu_files:
        print(f"No .vtu files found in {args.input_dir}")
        return

    print(f"Found {len(vtu_files)} frames in {args.input_dir}")

    # First pass: determine global value range for consistent colorbar.
    print("Scanning value range...")
    global_min = float("inf")
    global_max = float("-inf")
    for f in vtu_files:
        data = parse_vtu(f)
        if args.field not in data["cell_data"]:
            available = list(data["cell_data"].keys())
            print(f"Field '{args.field}' not found. Available: {available}")
            return
        vals = data["cell_data"][args.field]
        global_min = min(global_min, vals.min())
        global_max = max(global_max, vals.max())

    # Symmetric range for diverging colormap.
    vabs = max(abs(global_min), abs(global_max))
    if vabs == 0:
        vabs = 1.0
    print(f"Value range: [{-vabs:.4e}, {vabs:.4e}]")

    # Render frames.
    frames = []
    for i, f in enumerate(vtu_files):
        if i % 10 == 0:
            print(f"  rendering frame {i+1}/{len(vtu_files)}...")
        data = parse_vtu(f)
        img = render_frame(data, args.field, -vabs, vabs, i + 1, len(vtu_files))
        frames.append(img)

    # Save animation.
    output_path = args.output or str(Path(args.input_dir) / "animation.png")
    save_animation(frames, output_path, args.fps)
    fmt = output_format(output_path)
    if fmt == "gif":
        print("\nSaved {path} ({frames} frames, {fps} fps, GIF palette up to 256 colors)".format(
            path=output_path, frames=len(frames), fps=args.fps
        ))
    else:
        print("\nSaved {path} ({frames} frames, {fps} fps, full-color APNG)".format(
            path=output_path, frames=len(frames), fps=args.fps
        ))


if __name__ == "__main__":
    main()
