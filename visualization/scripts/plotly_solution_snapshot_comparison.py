#!/usr/bin/env python3
"""Plot deterministic and stochastic solution snapshots using Plotly."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PLOT_DIR = Path("visualization/plots")
DETERMINISTIC_SUFFIX = "_deterministic_path.csv"
STOCHASTIC_SUFFIX = "_stochastic_path.csv"
GOLDEN_RATIO = (1.0 + 5.0 ** 0.5) / 2.0
SCENE_DOMAIN_GAP = 0.012
RIGHT_COLORBAR_DOMAIN = 0.975
AXIS_TITLE_FONT_SIZE = 8
AXIS_TICK_FONT_SIZE = 7
COLORBAR_FONT_SIZE = 8
COLORBAR_TICK_FONT_SIZE = 7
DETERMINISTIC_Z_ASPECT_MULTIPLIER = 1.55


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an interactive Plotly comparison of the initial and final "
            "deterministic and stochastic solution snapshots."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Input snapshot CSV; defaults to the latest timestamped deterministic file.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output filename. If omitted, writes an HTML file under "
            "visualization/plots with a timestamped prefix."
        ),
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1280,
        help="Figure width in pixels.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=900,
        help="Figure height in pixels.",
    )
    parser.add_argument(
        "--colorscale",
        default="Cividis",
        help="Plotly colorscale for the surface panels.",
    )
    parser.add_argument(
        "--z-aspect",
        type=float,
        default=1.0 / GOLDEN_RATIO,
        help=(
            "Visual height of the value axis relative to the largest x/y axis. "
            "The default is 1/golden_ratio, preserving the square domain while "
            "making the xy plane visually dominant."
        ),
    )
    parser.add_argument(
        "--tracker",
        choices=["max-final-difference", "center", "none"],
        default="max-final-difference",
        help=(
            "Tracker marker location. 'max-final-difference' marks the grid "
            "point with largest |stochastic-deterministic| at final time."
        ),
    )
    return parser.parse_args()


def _derive_companion_path(given: Path) -> tuple[Path, Path]:
    name = given.name
    if name.endswith(STOCHASTIC_SUFFIX.lstrip("_")):
        det = given.parent / \
            name.replace("stochastic_path.csv", "deterministic_path.csv")
        return det, given
    if name.endswith(DETERMINISTIC_SUFFIX.lstrip("_")):
        sto = given.parent / \
            name.replace("deterministic_path.csv", "stochastic_path.csv")
        return given, sto
    raise ValueError(
        f"Input file must end with '{DETERMINISTIC_SUFFIX}' or '{STOCHASTIC_SUFFIX}': {given}"
    )


def resolve_input_paths(argument: str | None) -> tuple[Path, Path]:
    if argument is not None:
        return _derive_companion_path(Path(argument))

    candidates = sorted(DEFAULT_INPUT_DIR.glob(f"*{DETERMINISTIC_SUFFIX}"))
    if not candidates:
        raise FileNotFoundError(
            "No generated snapshot data found in data/output. First run 'make run-snapshot-comparison'."
        )
    det_path = candidates[-1]
    print(f"Using latest snapshot pair from {det_path}")
    return _derive_companion_path(det_path)


def _extract_timestamp_prefix(det_path: Path) -> str:
    stem = det_path.stem
    suffix = "_deterministic_path"
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def resolve_output_path(argument: str | None, det_path: Path) -> Path:
    if argument is not None:
        return Path(argument)
    prefix = _extract_timestamp_prefix(det_path)
    return DEFAULT_PLOT_DIR / f"{prefix}_solution_snapshot_comparison_plotly.html"


def _load_single_path_csv(
    path: Path,
) -> tuple[list[float], np.ndarray, np.ndarray, dict[float, np.ndarray], dict[str, str]]:
    expected_header = ["time", "x_index", "y_index", "x", "y", "value"]
    metadata: dict[str, str] = {}
    data_lines: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if stripped.startswith("#"):
                parts = stripped[1:].strip().split("=", maxsplit=1)
                if len(parts) == 2:
                    metadata[parts[0].strip()] = parts[1].strip()
                continue
            data_lines.append(raw_line)

    reader = csv.DictReader(data_lines)
    if reader.fieldnames != expected_header:
        raise ValueError(
            f"Expected CSV header {expected_header}, got {reader.fieldnames}")
    rows = list(reader)
    if not rows:
        raise ValueError(f"No snapshot rows found in {path}")

    times = sorted({float(row["time"]) for row in rows})
    nx = max(int(row["x_index"]) for row in rows)
    ny = max(int(row["y_index"]) for row in rows)
    x_coordinates = np.full(nx, np.nan, dtype=float)
    y_coordinates = np.full(ny, np.nan, dtype=float)
    fields: dict[float, np.ndarray] = {}
    for t in times:
        fields[t] = np.full((ny, nx), np.nan, dtype=float)

    for row in rows:
        t = float(row["time"])
        xi = int(row["x_index"]) - 1
        yi = int(row["y_index"]) - 1
        x_coordinates[xi] = float(row["x"])
        y_coordinates[yi] = float(row["y"])
        fields[t][yi, xi] = float(row["value"])

    if np.isnan(x_coordinates).any() or np.isnan(y_coordinates).any():
        raise ValueError(f"Incomplete coordinate coverage in {path}")
    for t, field in fields.items():
        if np.isnan(field).any():
            raise ValueError(
                f"Incomplete grid coverage for time={t} in {path}")

    return times, x_coordinates, y_coordinates, fields, metadata


def load_snapshot_data(
    det_path: Path,
    sto_path: Path,
) -> tuple[
    list[float],
    np.ndarray,
    np.ndarray,
    dict[tuple[str, float], np.ndarray],
    dict[str, str],
]:
    det_times, x_det, y_det, det_fields, metadata = _load_single_path_csv(
        det_path)
    sto_times, x_sto, y_sto, sto_fields, _ = _load_single_path_csv(sto_path)
    if det_times != sto_times:
        raise ValueError(
            "Deterministic and stochastic files have mismatched times.")
    if not np.array_equal(x_det, x_sto) or not np.array_equal(y_det, y_sto):
        raise ValueError(
            "Deterministic and stochastic files have mismatched coordinates.")

    fields: dict[tuple[str, float], np.ndarray] = {}
    for t in det_times:
        fields[("deterministic", t)] = det_fields[t]
        fields[("stochastic", t)] = sto_fields[t]
    return det_times, x_det, y_det, fields, metadata


def choose_color_scale(fields: dict[tuple[str, float], np.ndarray]) -> tuple[float, float, str]:
    all_values = np.concatenate([field.ravel() for field in fields.values()])
    vmin = float(np.min(all_values))
    vmax = float(np.max(all_values))
    if vmin < 0.0 < vmax:
        limit = max(abs(vmin), abs(vmax))
        return -limit, limit, "RdBu"
    return vmin, vmax, "Viridis"


def shared_value_range(fields: dict[tuple[str, float], np.ndarray]) -> tuple[float, float]:
    all_values = np.concatenate([field.ravel() for field in fields.values()])
    vmin = float(np.min(all_values))
    vmax = float(np.max(all_values))
    if np.isclose(vmin, vmax):
        padding = max(abs(vmin), 1.0) * 0.05
        return vmin - padding, vmax + padding
    padding = 0.03 * (vmax - vmin)
    return vmin - padding, vmax + padding


def fields_for_kind(
    display_fields: dict[tuple[str, float], np.ndarray],
    kind: str,
) -> dict[tuple[str, float], np.ndarray]:
    return {
        key: field
        for key, field in display_fields.items()
        if key[0] == kind
    }


def choose_tracker_index(
    display_fields: dict[tuple[str, float], np.ndarray],
    times: list[float],
    tracker_mode: str,
) -> tuple[int, int] | None:
    if tracker_mode == "none":
        return None

    sample = next(iter(display_fields.values()))
    ny, nx = sample.shape
    if tracker_mode == "center":
        return ny // 2, nx // 2

    final_time = times[-1]
    difference = (
        display_fields[("stochastic", final_time)]
        - display_fields[("deterministic", final_time)]
    )
    flat_index = int(np.nanargmax(np.abs(difference)))
    return tuple(int(index) for index in np.unravel_index(flat_index, difference.shape))


def format_time_value(time_value: float) -> str:
    return f"{time_value:g}"


def scene_domains() -> dict[int, dict[str, list[float]]]:
    x_mid = 0.5
    y_mid = 0.5
    half_gap = SCENE_DOMAIN_GAP / 2.0
    return {
        1: dict(
            x=[0.0, x_mid - half_gap],
            y=[y_mid + half_gap, 1.0],
        ),
        2: dict(
            x=[x_mid + half_gap, RIGHT_COLORBAR_DOMAIN],
            y=[y_mid + half_gap, 1.0],
        ),
        3: dict(
            x=[0.0, x_mid - half_gap],
            y=[0.0, y_mid - half_gap],
        ),
        4: dict(
            x=[x_mid + half_gap, RIGHT_COLORBAR_DOMAIN],
            y=[0.0, y_mid - half_gap],
        ),
    }


def make_plotly_figure(
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    display_fields: dict[tuple[str, float], np.ndarray],
    times: list[float],
    colorscale: str,
    width: int,
    height: int,
    z_aspect: float,
    tracker_mode: str,
) -> go.Figure:
    figure = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"type": "surface"}, {"type": "surface"}],
               [{"type": "surface"}, {"type": "surface"}]],
        vertical_spacing=SCENE_DOMAIN_GAP,
        horizontal_spacing=SCENE_DOMAIN_GAP,
    )

    panel_order = [
        ("deterministic", times[0], 1, 1),
        ("stochastic", times[0], 1, 2),
        ("deterministic", times[-1], 2, 1),
        ("stochastic", times[-1], 2, 2),
    ]

    kind_color_scales = {
        kind: choose_color_scale(fields_for_kind(display_fields, kind))
        for kind in ["deterministic", "stochastic"]
    }
    kind_z_ranges = {
        kind: shared_value_range(fields_for_kind(display_fields, kind))
        for kind in ["deterministic", "stochastic"]
    }

    x_range = float(np.nanmax(x_coordinates) - np.nanmin(x_coordinates))
    y_range = float(np.nanmax(y_coordinates) - np.nanmin(y_coordinates))
    floor_max = max(x_range, y_range, 1.0)
    aspectratio = dict(x=x_range / floor_max,
                       y=y_range / floor_max, z=z_aspect)
    deterministic_aspectratio = dict(
        x=aspectratio["x"],
        y=aspectratio["y"],
        z=z_aspect * DETERMINISTIC_Z_ASPECT_MULTIPLIER,
    )
    tracker_index = choose_tracker_index(display_fields, times, tracker_mode)
    surface_trace_indices: list[int] = []

    for kind, time_value, row, col in panel_order:
        field = display_fields[(kind, time_value)]
        vmin, vmax, default_colorscale = kind_color_scales[kind]
        colorscale_name = colorscale or default_colorscale
        colorbar_x = 0.49 if kind == "deterministic" else 0.99
        trace = go.Surface(
            x=x_coordinates,
            y=y_coordinates,
            z=field,
            surfacecolor=field,
            colorscale=colorscale_name,
            cmin=vmin,
            cmax=vmax,
            showscale=(row == 2),
            opacity=0.94,
            colorbar=dict(
                title=dict(
                    text=f"{kind.capitalize()} scale",
                    font=dict(size=COLORBAR_FONT_SIZE),
                ),
                lenmode="fraction",
                len=0.42,
                x=colorbar_x,
                y=0.24,
                thickness=14,
                tickfont=dict(size=COLORBAR_TICK_FONT_SIZE),
            ),
            contours_z=dict(
                show=True,
                usecolormap=True,
                highlightcolor="black",
                highlightwidth=1,
                project_z=True,
            ),
            hovertemplate=(
                f"{kind.capitalize()}<br>"
                f"t={format_time_value(time_value)}<br>"
                "x=%{x:.4g}<br>y=%{y:.4g}<br>"
                "value=%{z:.6g}<extra></extra>"
            ),
        )
        figure.add_trace(trace, row=row, col=col)
        surface_trace_indices.append(len(figure.data) - 1)

    if tracker_index is not None:
        tracker_y_index, tracker_x_index = tracker_index
        tracker_x = float(x_coordinates[tracker_x_index])
        tracker_y = float(y_coordinates[tracker_y_index])
        for kind, time_value, row, col in panel_order:
            field = display_fields[(kind, time_value)]
            deterministic_value = float(
                display_fields[("deterministic", time_value)][
                    tracker_y_index, tracker_x_index
                ]
            )
            stochastic_value = float(
                display_fields[("stochastic", time_value)][
                    tracker_y_index, tracker_x_index
                ]
            )
            difference = stochastic_value - deterministic_value
            z_value = float(field[tracker_y_index, tracker_x_index])
            marker = go.Scatter3d(
                x=[tracker_x],
                y=[tracker_y],
                z=[z_value],
                mode="markers",
                marker=dict(
                    size=5,
                    color="black",
                    line=dict(color="yellow", width=3),
                    symbol="diamond",
                ),
                name="tracker",
                showlegend=False,
                hovertemplate=(
                    f"Tracker ({kind})<br>"
                    f"t={format_time_value(time_value)}<br>"
                    f"x={tracker_x:.6g}<br>"
                    f"y={tracker_y:.6g}<br>"
                    f"deterministic={deterministic_value:.6g}<br>"
                    f"stochastic={stochastic_value:.6g}<br>"
                    f"stochastic-deterministic={difference:.6g}"
                    "<extra></extra>"
                ),
            )
            figure.add_trace(marker, row=row, col=col)

    camera_golden = dict(eye=dict(x=1.32, y=-1.48, z=0.6))
    camera_iso = dict(eye=dict(x=1.25, y=-1.25, z=0.86))
    camera_top = dict(eye=dict(x=0.0, y=0.0, z=3.0))
    camera_side = dict(eye=dict(x=2.2, y=0.0, z=0.36))
    domains = scene_domains()

    scene_template = dict(
        xaxis=dict(
            title=dict(text="x", font=dict(size=AXIS_TITLE_FONT_SIZE)),
            showgrid=True,
            zeroline=False,
            range=[float(np.nanmin(x_coordinates)), float(np.nanmax(x_coordinates))],
            tickfont=dict(size=AXIS_TICK_FONT_SIZE),
            nticks=4,
        ),
        yaxis=dict(
            title=dict(text="y", font=dict(size=AXIS_TITLE_FONT_SIZE)),
            showgrid=True,
            zeroline=False,
            range=[float(np.nanmin(y_coordinates)), float(np.nanmax(y_coordinates))],
            tickfont=dict(size=AXIS_TICK_FONT_SIZE),
            nticks=4,
        ),
        camera=camera_golden,
        aspectmode="manual",
        bgcolor="rgba(0,0,0,0)",
    )

    scene_kinds = {
        1: "deterministic",
        2: "stochastic",
        3: "deterministic",
        4: "stochastic",
    }
    for scene_id, scene_kind in scene_kinds.items():
        zmin, zmax = kind_z_ranges[scene_kind]
        scene_layout = dict(scene_template)
        scene_layout["domain"] = domains[scene_id]
        scene_layout["aspectratio"] = (
            deterministic_aspectratio
            if scene_kind == "deterministic"
            else aspectratio
        )
        scene_layout["zaxis"] = dict(
            title=dict(
                text=(
                    "value (visual z x1.55)"
                    if scene_kind == "deterministic"
                    else "value"
                ),
                font=dict(size=AXIS_TITLE_FONT_SIZE),
            ),
            showgrid=True,
            zeroline=False,
            range=[zmin, zmax],
            tickfont=dict(size=AXIS_TICK_FONT_SIZE),
            nticks=5,
        )
        figure.update_layout({f"scene{scene_id}": scene_layout})

    colorscale_buttons = [
        dict(label=scale, method="restyle", args=[{"colorscale": scale}])
        for scale in ["Cividis", "Viridis", "Plasma", "RdBu"]
    ]

    view_buttons = [
        dict(
            label="Gold",
            method="relayout",
            args=[
                {
                    "scene.camera": camera_golden,
                    "scene2.camera": camera_golden,
                    "scene3.camera": camera_golden,
                    "scene4.camera": camera_golden,
                }
            ],
        ),
        dict(
            label="Iso",
            method="relayout",
            args=[
                {
                    "scene.camera": camera_iso,
                    "scene2.camera": camera_iso,
                    "scene3.camera": camera_iso,
                    "scene4.camera": camera_iso,
                }
            ],
        ),
        dict(
            label="↑",
            method="relayout",
            args=[
                {
                    "scene.camera": camera_top,
                    "scene2.camera": camera_top,
                    "scene3.camera": camera_top,
                    "scene4.camera": camera_top,
                }
            ],
        ),
        dict(
            label="→",
            method="relayout",
            args=[
                {
                    "scene.camera": camera_side,
                    "scene2.camera": camera_side,
                    "scene3.camera": camera_side,
                    "scene4.camera": camera_side,
                }
            ],
        ),
    ]

    contour_buttons = [
        dict(
            label="On",
            method="restyle",
            args=[{"contours.z.show": True}, surface_trace_indices],
        ),
        dict(
            label="Off",
            method="restyle",
            args=[{"contours.z.show": False}, surface_trace_indices],
        ),
    ]

    figure.update_layout(
        title=dict(
            text="Initial/Final Deterministic/Stochastic Snapshots",
            font=dict(size=14),
            x=0.5,
            xanchor="center",
        ),
        width=width,
        height=height,
        margin=dict(l=0, r=48, t=44, b=0),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                buttons=colorscale_buttons,
                x=0.01,
                y=1.045,
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1,
                font=dict(size=10),
            ),
            dict(
                type="buttons",
                direction="right",
                buttons=view_buttons,
                x=0.38,
                y=1.045,
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1,
                font=dict(size=10),
            ),
            dict(
                type="buttons",
                direction="right",
                buttons=contour_buttons,
                x=0.68,
                y=1.045,
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1,
                font=dict(size=10),
            ),
        ],
    )

    figure.add_annotation(
        x=0.19,
        y=1.015,
        xref="paper",
        yref="paper",
        text="<b>Colorscale</b>",
        showarrow=False,
        font=dict(size=9),
    )
    figure.add_annotation(
        x=0.56,
        y=1.015,
        xref="paper",
        yref="paper",
        text="<b>Camera</b>",
        showarrow=False,
        font=dict(size=9),
    )
    figure.add_annotation(
        x=0.85,
        y=1.015,
        xref="paper",
        yref="paper",
        text="<b>Contours</b>",
        showarrow=False,
        font=dict(size=9),
    )
    if tracker_index is not None:
        tracker_y_index, tracker_x_index = tracker_index
        tracker_x = float(x_coordinates[tracker_x_index])
        tracker_y = float(y_coordinates[tracker_y_index])
        final_difference = float(
            display_fields[("stochastic", times[-1])][tracker_y_index, tracker_x_index]
            - display_fields[("deterministic", times[-1])][
                tracker_y_index, tracker_x_index
            ]
        )
        figure.add_annotation(
            x=0.99,
            y=1.005,
            xref="paper",
            yref="paper",
            xanchor="right",
            yanchor="bottom",
            text=(
                f"Max noise effect: (x={tracker_x:.3g}, y={tracker_y:.3g}), "
                f"Δ={final_difference:.3g}"
            ),
            showarrow=False,
            font=dict(size=8),
            bgcolor="rgba(255,255,255,0.6)",
            borderpad=2,
        )

    return figure


def main() -> None:
    args = parse_args()
    det_path, sto_path = resolve_input_paths(args.input)
    output_path = resolve_output_path(args.output, det_path)

    times, x_coordinates, y_coordinates, fields, _ = load_snapshot_data(
        det_path, sto_path)
    display_times = [times[0], times[-1]]
    display_fields = {
        (kind, t): fields[(kind, t)]
        for kind in ["deterministic", "stochastic"]
        for t in display_times
    }

    figure = make_plotly_figure(
        x_coordinates,
        y_coordinates,
        display_fields,
        display_times,
        args.colorscale,
        args.width,
        args.height,
        args.z_aspect,
        args.tracker,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".html":
        figure.write_html(
            str(output_path),
            include_plotlyjs="cdn",
            config={
                "responsive": True,
                "displayModeBar": True,
                "modeBarButtonsToRemove": ["select2d", "lasso2d"],
            },
        )
    else:
        figure.write_image(str(output_path))

    print(f"Wrote Plotly snapshot comparison to {output_path}")


if __name__ == "__main__":
    main()
