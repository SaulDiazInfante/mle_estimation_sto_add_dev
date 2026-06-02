#!/usr/bin/env python3
"""Plot reconstructed deterministic and stochastic solution snapshots."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PLOT_DIR = Path("visualization/plots")
SNAPSHOT_SUFFIX = "_solution_snapshot_comparison.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot deterministic and stochastic solution snapshots at matching times."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Input snapshot CSV; defaults to the latest timestamped snapshot file.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output image path; defaults to a matching timestamped PNG.",
    )
    parser.add_argument(
        "--velocity-mode-x",
        type=int,
        default=None,
        help="Override the x-mode of the stationary velocity field.",
    )
    parser.add_argument(
        "--velocity-mode-y",
        type=int,
        default=None,
        help="Override the y-mode of the stationary velocity field.",
    )
    parser.add_argument(
        "--length-x",
        type=float,
        default=None,
        help="Override the x-domain length used in the velocity field.",
    )
    parser.add_argument(
        "--length-y",
        type=float,
        default=None,
        help="Override the y-domain length used in the velocity field.",
    )
    parser.add_argument(
        "--quiver-stride",
        type=int,
        default=None,
        help="Plot every k-th velocity vector in each direction.",
    )
    return parser.parse_args()


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


def resolve_output_path(argument: str | None, input_path: Path) -> Path:
    if argument is not None:
        return Path(argument)

    output_name = f"{input_path.stem}.png" if input_path.suffix == ".csv" else f"{input_path.name}.png"
    return DEFAULT_PLOT_DIR / output_name


def load_snapshot_data(
    path: Path,
) -> tuple[
    list[float],
    np.ndarray,
    np.ndarray,
    dict[tuple[str, float], np.ndarray],
    dict[str, str],
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
    metadata: dict[str, str] = {}
    rows: list[dict[str, str]] = []

    with path.open("r", encoding="utf-8") as handle:
        data_lines: list[str] = []
        for raw_line in handle:
            stripped_line = raw_line.strip()
            if stripped_line.startswith("#"):
                key_value = stripped_line[1:].strip().split("=", maxsplit=1)
                if len(key_value) == 2:
                    metadata[key_value[0].strip()] = key_value[1].strip()
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

    return times, x_coordinates, y_coordinates, fields, metadata


def centers_to_edges(centers: np.ndarray) -> np.ndarray:
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=float)

    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def format_time_label(time_value: float) -> str:
    return f"{time_value:g}"


def choose_color_scale(fields: dict[tuple[str, float], np.ndarray]) -> tuple[float, float, str]:
    all_values = np.concatenate([field.ravel() for field in fields.values()])
    value_min = float(np.min(all_values))
    value_max = float(np.max(all_values))

    if value_min < 0.0 < value_max:
        limit = max(abs(value_min), abs(value_max))
        return -limit, limit, "RdBu_r"

    return value_min, value_max, "viridis"


def read_int_setting(
    explicit_value: int | None,
    metadata: dict[str, str],
    metadata_key: str,
    env_name: str,
    default_value: int,
) -> int:
    if explicit_value is not None:
        return explicit_value
    if metadata_key in metadata:
        return int(metadata[metadata_key])
    if env_name in os.environ:
        return int(os.environ[env_name])
    return default_value


def read_float_setting(
    explicit_value: float | None,
    metadata: dict[str, str],
    metadata_key: str,
    env_name: str,
    default_value: float,
) -> float:
    if explicit_value is not None:
        return explicit_value
    if metadata_key in metadata:
        return float(metadata[metadata_key])
    if env_name in os.environ:
        return float(os.environ[env_name])
    return default_value


def compute_velocity_field(
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    velocity_mode_x: int,
    velocity_mode_y: int,
    length_x: float,
    length_y: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_grid, y_grid = np.meshgrid(x_coordinates, y_coordinates)
    wave_x = np.pi * velocity_mode_x * x_grid / length_x
    wave_y = np.pi * velocity_mode_y * y_grid / length_y
    velocity_x = (
        -2.0
        * np.pi
        * velocity_mode_y
        / length_y
        * np.sin(wave_x)
        * np.cos(wave_y)
    )
    velocity_y = (
        2.0
        * np.pi
        * velocity_mode_x
        / length_x
        * np.cos(wave_x)
        * np.sin(wave_y)
    )
    return x_grid, y_grid, velocity_x, velocity_y


def add_velocity_quiver(
    ax: plt.Axes,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    velocity_x: np.ndarray,
    velocity_y: np.ndarray,
    stride: int,
) -> None:
    ax.quiver(
        x_grid[::stride, ::stride],
        y_grid[::stride, ::stride],
        velocity_x[::stride, ::stride],
        velocity_y[::stride, ::stride],
        angles="xy",
        scale_units="xy",
        scale=None,
        pivot="mid",
        color="#111111",
        alpha=0.45,
        width=0.0045,
        headwidth=3.4,
        headlength=4.6,
        headaxislength=4.0,
        zorder=3,
    )


def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)
    output_path = resolve_output_path(args.output, input_path)

    times, x_coordinates, y_coordinates, fields, metadata = load_snapshot_data(input_path)
    x_edges = centers_to_edges(x_coordinates)
    y_edges = centers_to_edges(y_coordinates)
    vmin, vmax, cmap = choose_color_scale(fields)
    velocity_mode_x = read_int_setting(
        args.velocity_mode_x,
        metadata,
        "velocity_mode_x",
        "SPDE_VELOCITY_MODE_X",
        1,
    )
    velocity_mode_y = read_int_setting(
        args.velocity_mode_y,
        metadata,
        "velocity_mode_y",
        "SPDE_VELOCITY_MODE_Y",
        1,
    )
    length_x = read_float_setting(
        args.length_x,
        metadata,
        "length_x",
        "SPDE_LENGTH_X",
        5.0,
    )
    length_y = read_float_setting(
        args.length_y,
        metadata,
        "length_y",
        "SPDE_LENGTH_Y",
        5.0,
    )
    quiver_stride = args.quiver_stride
    if quiver_stride is None:
        quiver_stride = max(1, int(np.ceil(max(len(x_coordinates), len(y_coordinates)) / 16.0)))
    x_grid, y_grid, velocity_x, velocity_y = compute_velocity_field(
        x_coordinates,
        y_coordinates,
        velocity_mode_x,
        velocity_mode_y,
        length_x,
        length_y,
    )

    fig, axes = plt.subplots(
        2,
        len(times),
        figsize=(4.1 * len(times), 6.3),
        sharex=True,
        sharey=True,
    )
    if len(times) == 1:
        axes = np.array(axes).reshape(2, 1)

    mesh = None
    row_labels = ["Deterministic", "Stochastic"]
    row_kinds = ["deterministic", "stochastic"]
    for row_index, kind in enumerate(row_kinds):
        for column_index, time_value in enumerate(times):
            ax = axes[row_index, column_index]
            mesh = ax.pcolormesh(
                x_edges,
                y_edges,
                fields[(kind, time_value)],
                shading="auto",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            add_velocity_quiver(
                ax,
                x_grid,
                y_grid,
                velocity_x,
                velocity_y,
                quiver_stride,
            )
            ax.set_aspect("equal")
            ax.set_title(rf"$t = {format_time_label(time_value)}$")
            if row_index == len(row_kinds) - 1:
                ax.set_xlabel(r"$x$")
            if column_index == 0:
                ax.set_ylabel(r"$y$")
            ax.tick_params(direction="out", length=3.0)

    fig.text(0.03, 0.70, row_labels[0], rotation=90, va="center", ha="center", fontsize=12)
    fig.text(0.03, 0.28, row_labels[1], rotation=90, va="center", ha="center", fontsize=12)
    colorbar = fig.colorbar(mesh, ax=axes, fraction=0.030, pad=0.040)
    colorbar.set_label("Solution value")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(
        left=0.08,
        right=0.86,
        bottom=0.10,
        top=0.92,
        wspace=0.14,
        hspace=0.18,
    )
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
