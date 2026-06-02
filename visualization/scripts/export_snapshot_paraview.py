#!/usr/bin/env python3
"""Export reconstructed solution snapshots to ParaView VTK XML files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PARAVIEW_DIR = Path("visualization/paraview")
DEFAULT_COLLECTION_NAME = "solution_snapshots"
SNAPSHOT_SUFFIX = "_solution_snapshot_comparison.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export deterministic and stochastic solution snapshots as a "
            "ParaView .pvd time collection with .vts structured-grid frames."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Input snapshot CSV; defaults to the latest timestamped snapshot file.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Defaults to visualization/paraview/<input-stem>_paraview."
        ),
    )
    parser.add_argument(
        "--collection-name",
        default=DEFAULT_COLLECTION_NAME,
        help="Base name for the generated .pvd collection.",
    )
    parser.add_argument(
        "--warp-scale",
        type=float,
        default=0.0,
        help=(
            "Scale factor applied to the stochastic field to set z-coordinates. "
            "0 (default) keeps the flat z=0 plane. "
            "A non-zero value lifts the surface into 3-D so ParaView can rotate it."
        ),
    )
    return parser.parse_args()


def resolve_output_dir(argument: str | None, input_path: Path) -> Path:
    if argument is not None:
        return Path(argument)

    stem = input_path.stem
    suffix = "_solution_snapshot_comparison"
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    return DEFAULT_PARAVIEW_DIR / f"{stem}_paraview"


def resolve_input_path(argument: str | None) -> Path:
    if argument is not None:
        return Path(argument)

    candidates = sorted(DEFAULT_INPUT_DIR.glob(f"*{SNAPSHOT_SUFFIX}"))
    if not candidates:
        msg = (
            "No generated snapshot data was found in "
            f"{DEFAULT_INPUT_DIR}. First run 'make run-snapshot-comparison'."
        )
        raise FileNotFoundError(msg)

    print(f"Using latest snapshot comparison from {candidates[-1]}")
    return candidates[-1]


def load_snapshot_data(
    path: Path,
) -> tuple[
    list[float],
    np.ndarray,
    np.ndarray,
    dict[tuple[str, float], np.ndarray],
]:
    expected_header = [
        "solution_kind",
        "time",
        "x_index",
        "y_index",
        "x",
        "y",
        "value",
    ]
    rows: list[dict[str, str]] = []

    with path.open("r", encoding="utf-8") as handle:
        data_lines: list[str] = []
        for raw_line in handle:
            if raw_line.strip().startswith("#"):
                continue
            data_lines.append(raw_line)

        reader = csv.DictReader(data_lines)
        if reader.fieldnames != expected_header:
            msg = f"Expected CSV header {expected_header}, got {reader.fieldnames}"
            raise ValueError(msg)
        rows = list(reader)

    if not rows:
        msg = f"No snapshot rows found in {path}"
        raise ValueError(msg)

    times = sorted({float(row["time"]) for row in rows})
    kinds = sorted({row["solution_kind"] for row in rows})
    if kinds != ["deterministic", "stochastic"]:
        msg = f"Expected solution kinds ['deterministic', 'stochastic'], got {kinds}"
        raise ValueError(msg)

    nx = max(int(row["x_index"]) for row in rows)
    ny = max(int(row["y_index"]) for row in rows)
    x_coordinates = np.full(nx, np.nan, dtype=float)
    y_coordinates = np.full(ny, np.nan, dtype=float)

    fields: dict[tuple[str, float], np.ndarray] = {}
    for kind in kinds:
        for time_value in times:
            fields[(kind, time_value)] = np.full((ny, nx), np.nan, dtype=float)

    for row in rows:
        kind = row["solution_kind"]
        time_value = float(row["time"])
        x_index = int(row["x_index"]) - 1
        y_index = int(row["y_index"]) - 1
        x_coordinates[x_index] = float(row["x"])
        y_coordinates[y_index] = float(row["y"])
        fields[(kind, time_value)][y_index, x_index] = float(row["value"])

    if np.isnan(x_coordinates).any() or np.isnan(y_coordinates).any():
        msg = f"Incomplete coordinate coverage in {path}"
        raise ValueError(msg)

    for key, field in fields.items():
        if np.isnan(field).any():
            msg = f"Incomplete grid coverage for {key} in {path}"
            raise ValueError(msg)

    return times, x_coordinates, y_coordinates, fields


def write_data_array(
    handle,
    name: str,
    values: np.ndarray,
    components: int = 1,
) -> None:
    component_text = ""
    if components > 1:
        component_text = f' NumberOfComponents="{components}"'
    handle.write(
        f'        <DataArray type="Float64" Name="{escape(name)}"'
        f'{component_text} format="ascii">\n'
    )

    flat_values = values.reshape(-1, components)
    for row in flat_values:
        handle.write("          ")
        handle.write(" ".join(f"{float(value):.17e}" for value in row))
        handle.write("\n")

    handle.write("        </DataArray>\n")


def write_vts_frame(
    frame_path: Path,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    deterministic_field: np.ndarray,
    stochastic_field: np.ndarray,
    warp_scale: float = 0.0,
) -> None:
    nx = len(x_coordinates)
    ny = len(y_coordinates)
    x_grid, y_grid = np.meshgrid(x_coordinates, y_coordinates)

    if warp_scale != 0.0:
        # Lift the surface into 3-D using the stochastic field as z.
        # This avoids ParaView's 2-D interaction lock on flat single-layer grids.
        z_grid = warp_scale * stochastic_field
        nz = len(np.unique(z_grid))
        extent = f"0 {nx - 1} 0 {ny - 1} 0 {max(nz - 1, 1)}"
    else:
        z_grid = np.zeros_like(x_grid)
        extent = f"0 {nx - 1} 0 {ny - 1} 0 0"
    points = np.stack([x_grid, y_grid, z_grid], axis=-1)

    with frame_path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write(
            '<VTKFile type="StructuredGrid" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        handle.write(f'  <StructuredGrid WholeExtent="{extent}">\n')
        handle.write(f'    <Piece Extent="{extent}">\n')
        handle.write('      <PointData Scalars="stochastic">\n')
        write_data_array(handle, "deterministic", deterministic_field)
        write_data_array(handle, "stochastic", stochastic_field)
        handle.write("      </PointData>\n")
        handle.write("      <CellData/>\n")
        handle.write("      <Points>\n")
        write_data_array(handle, "points", points, components=3)
        handle.write("      </Points>\n")
        handle.write("    </Piece>\n")
        handle.write("  </StructuredGrid>\n")
        handle.write("</VTKFile>\n")


def write_pvd_collection(
    collection_path: Path,
    frame_paths: list[Path],
    times: list[float],
) -> None:
    with collection_path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write(
            '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        )
        handle.write("  <Collection>\n")
        for time_value, frame_path in zip(times, frame_paths):
            relative_path = frame_path.relative_to(collection_path.parent)
            handle.write(
                f'    <DataSet timestep="{time_value:.17e}" group="" part="0" '
                f'file="{escape(relative_path.as_posix())}"/>\n'
            )
        handle.write("  </Collection>\n")
        handle.write("</VTKFile>\n")


def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)
    output_dir = resolve_output_dir(args.output_dir, input_path)
    frame_dir = output_dir / "frames"

    times, x_coordinates, y_coordinates, fields = load_snapshot_data(input_path)

    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[Path] = []

    for frame_index, time_value in enumerate(times):
        frame_path = frame_dir / f"solution_{frame_index:04d}.vts"
        write_vts_frame(
            frame_path,
            x_coordinates,
            y_coordinates,
            fields[("deterministic", time_value)],
            fields[("stochastic", time_value)],
            warp_scale=args.warp_scale,
        )
        frame_paths.append(frame_path)

    collection_path = output_dir / f"{args.collection_name}.pvd"
    write_pvd_collection(collection_path, frame_paths, times)

    print(f"Wrote ParaView collection to {collection_path}")
    print(f"Wrote {len(frame_paths)} VTK structured-grid frame(s) to {frame_dir}")


if __name__ == "__main__":
    main()
