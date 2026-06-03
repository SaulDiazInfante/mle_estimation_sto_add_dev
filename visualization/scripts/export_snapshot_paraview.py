#!/usr/bin/env python3
"""Export reconstructed solution snapshots to ParaView VTK XML files."""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PARAVIEW_DIR = Path("visualization/paraview")
DEFAULT_COLLECTION_NAME = "solution_snapshots"
DETERMINISTIC_SUFFIX = "_deterministic_path.csv"
STOCHASTIC_SUFFIX = "_stochastic_path.csv"


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


def _derive_companion_path(given: Path) -> tuple[Path, Path]:
    """Return (deterministic_path, stochastic_path) given either file."""
    name = given.name
    if name.endswith(STOCHASTIC_SUFFIX.lstrip("_")):
        det = given.parent / name.replace("stochastic_path.csv", "deterministic_path.csv")
        return det, given
    if name.endswith(DETERMINISTIC_SUFFIX.lstrip("_")):
        sto = given.parent / name.replace("deterministic_path.csv", "stochastic_path.csv")
        return given, sto
    msg = (
        f"Input file must end with '{DETERMINISTIC_SUFFIX}' or "
        f"'{STOCHASTIC_SUFFIX}': {given}"
    )
    raise ValueError(msg)


def _extract_timestamp_prefix(det_path: Path) -> str:
    stem = det_path.stem
    suffix = "_deterministic_path"
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def resolve_output_dir(argument: str | None, det_path: Path) -> Path:
    if argument is not None:
        return Path(argument)

    prefix = _extract_timestamp_prefix(det_path)
    return DEFAULT_PARAVIEW_DIR / f"{prefix}_paraview"


def resolve_input_paths(argument: str | None) -> tuple[Path, Path]:
    if argument is not None:
        return _derive_companion_path(Path(argument))

    candidates = sorted(DEFAULT_INPUT_DIR.glob(f"*{DETERMINISTIC_SUFFIX}"))
    if not candidates:
        msg = (
            "No generated snapshot data was found in "
            f"{DEFAULT_INPUT_DIR}. First run 'make run-snapshot-comparison'."
        )
        raise FileNotFoundError(msg)

    det_path = candidates[-1]
    print(f"Using latest snapshot pair from {det_path}")
    return _derive_companion_path(det_path)


def _stream_single_path(path: Path):
    """Yield (time, x_coordinates, y_coordinates, field) per snapshot from one file."""
    expected_header = ["time", "x_index", "y_index", "x", "y", "value"]

    def non_comment_lines(handle):
        for raw_line in handle:
            if not raw_line.strip().startswith("#"):
                yield raw_line

    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(non_comment_lines(handle))
        if list(reader.fieldnames) != expected_header:
            msg = f"Expected CSV header {expected_header}, got {reader.fieldnames}"
            raise ValueError(msg)

        for time_key, group_iter in itertools.groupby(
            reader, key=lambda r: r["time"]
        ):
            time_value = float(time_key)
            group = list(group_iter)

            nx = max(int(r["x_index"]) for r in group)
            ny = max(int(r["y_index"]) for r in group)
            x_coordinates = np.full(nx, np.nan, dtype=float)
            y_coordinates = np.full(ny, np.nan, dtype=float)
            field = np.full((ny, nx), np.nan, dtype=float)

            for row in group:
                xi = int(row["x_index"]) - 1
                yi = int(row["y_index"]) - 1
                x_coordinates[xi] = float(row["x"])
                y_coordinates[yi] = float(row["y"])
                field[yi, xi] = float(row["value"])

            yield time_value, x_coordinates, y_coordinates, field


def stream_snapshot_frames(det_path: Path, sto_path: Path):
    """Yield (time, x_coordinates, y_coordinates, det_field, sto_field)
    by reading both single-kind CSVs in lockstep."""
    for det_frame, sto_frame in zip(
        _stream_single_path(det_path), _stream_single_path(sto_path)
    ):
        time_value, x_coords, y_coords, det_field = det_frame
        _, _, _, sto_field = sto_frame
        yield time_value, x_coords, y_coords, det_field, sto_field


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

    # The structured grid always has a single z-layer (extent 0 0).
    # When warp_scale != 0, z-coordinates carry the field values so the
    # surface has depth in ParaView without changing the grid topology.
    extent = f"0 {nx - 1} 0 {ny - 1} 0 0"
    if warp_scale != 0.0:
        z_grid = warp_scale * stochastic_field
    else:
        z_grid = np.zeros_like(x_grid)
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
    det_path, sto_path = resolve_input_paths(args.input)
    output_dir = resolve_output_dir(args.output_dir, det_path)
    frame_dir = output_dir / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)

    frame_paths: list[Path] = []
    times: list[float] = []
    report_every = 1000  # print progress every N frames

    print(f"Streaming frames from {det_path} and {sto_path} ...")
    for frame_index, (time_value, x_coords, y_coords, det_field, sto_field) in enumerate(
        stream_snapshot_frames(det_path, sto_path)
    ):
        frame_path = frame_dir / f"solution_{frame_index:05d}.vts"
        write_vts_frame(
            frame_path,
            x_coords,
            y_coords,
            det_field,
            sto_field,
            warp_scale=args.warp_scale,
        )
        frame_paths.append(frame_path)
        times.append(time_value)

        if (frame_index + 1) % report_every == 0:
            print(f"  {frame_index + 1} frames written ...", flush=True)

    collection_path = output_dir / f"{args.collection_name}.pvd"
    write_pvd_collection(collection_path, frame_paths, times)

    print(f"Wrote ParaView collection to {collection_path}")
    print(f"Wrote {len(frame_paths)} VTK structured-grid frame(s) to {frame_dir}")


if __name__ == "__main__":
    main()
