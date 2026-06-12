#!/usr/bin/env python3
"""Export manuscript-ready 3-D snapshot panels."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS_DIR))
from plot_solution_snapshot_comparison import (  # noqa: E402
    DEFAULT_PLOT_DIR,
    _extract_timestamp_prefix,
    resolve_input_paths,
)


EXPECTED_SNAPSHOT_HEADER = ["time", "x_index", "y_index", "x", "y", "value"]


PANEL_SPECS = [
    ("initial", "deterministic", 0),
    ("final", "deterministic", -1),
    ("initial", "stochastic", 0),
    ("final", "stochastic", -1),
]

PANEL_IDS = ["(a)", "(b)", "(c)", "(d)"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write high-resolution 3-D surface panels for the initial/final "
            "deterministic/stochastic snapshots."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Snapshot CSV; defaults to the latest timestamped file.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for exported figures. Defaults to "
            "visualization/plots/<timestamp>_manuscript_panels."
        ),
    )
    parser.add_argument(
        "--extension",
        default="png",
        choices=["png", "tif", "tiff", "jpg", "jpeg"],
        help="Image format. PNG is convenient; TIFF is recommended for submission.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=500,
        help="Output resolution in dots per inch.",
    )
    parser.add_argument(
        "--width-inches",
        type=float,
        default=3.3,
        help="Single-panel figure width (Elsevier 1-column standard is 3.3 in).",
    )
    parser.add_argument(
        "--height-inches",
        type=float,
        default=3.0,
        help="Single-panel figure height.",
    )
    parser.add_argument(
        "--elev",
        type=float,
        default=29.0,
        help="Camera elevation angle in degrees.",
    )
    parser.add_argument(
        "--azim",
        type=float,
        default=225.0,
        help="Camera azimuth angle in degrees.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Mesh stride for the rendered surface.",
    )
    parser.add_argument(
        "--contour-levels",
        type=int,
        default=12,
        help="Number of floor contour levels.",
    )
    parser.add_argument(
        "--cmap",
        default=None,
        help="Global colormap override. If None, uses kind-specific defaults.",
    )
    parser.add_argument(
        "--cmap-deterministic",
        default="turbo",
        help="Colormap for deterministic panels.",
    )
    parser.add_argument(
        "--cmap-stochastic",
        default="Spectral_r",
        help="Colormap for stochastic panels.",
    )
    parser.add_argument(
        "--deterministic-z-boost",
        type=float,
        default=1.15,
        help="Visual z-axis boost for deterministic panels.",
    )
    parser.add_argument(
        "--stochastic-z-boost",
        type=float,
        default=1.15,
        help="Visual z-axis boost for stochastic panels.",
    )
    parser.add_argument(
        "--variant",
        choices=["signed", "nonnegative", "both"],
        default="both",
        help="Export signed values, negative-masked values, or both.",
    )
    parser.add_argument(
        "--layout",
        choices=["separate", "combined", "both"],
        default="both",
        help="Export separate panel files, one combined 2x2 figure, or both.",
    )
    parser.add_argument(
        "--enhanced-floor",
        action="store_true",
        help="Increase line thickness of floor projections.",
    )
    parser.add_argument(
        "--deterministic-inset",
        action="store_true",
        help=(
            "Deprecated compatibility flag. Insets are now enabled for every "
            "panel by default; use --no-insets to suppress them."
        ),
    )
    parser.add_argument(
        "--insets",
        dest="insets",
        action="store_true",
        default=True,
        help="Add local zoom insets centered on each panel's largest magnitude.",
    )
    parser.add_argument(
        "--no-insets",
        dest="insets",
        action="store_false",
        help="Do not draw local zoom insets.",
    )
    parser.add_argument(
        "--inset-zoom-fraction",
        type=float,
        default=0.34,
        help="Fraction of each spatial axis shown in the local zoom inset.",
    )
    parser.add_argument(
        "--surface-wireframe",
        dest="surface_wireframe",
        action="store_true",
        default=True,
        help="Overlay a subtle black mesh on surfaces to clarify low-relief panels.",
    )
    parser.add_argument(
        "--no-surface-wireframe",
        dest="surface_wireframe",
        action="store_false",
        help="Do not overlay a black mesh on surfaces.",
    )
    parser.add_argument(
        "--wireframe-stride",
        type=int,
        default=4,
        help="Mesh stride for the black surface wireframe overlay.",
    )
    parser.add_argument(
        "--wireframe-linewidth",
        type=float,
        default=0.18,
        help="Line width for the black surface wireframe overlay.",
    )
    parser.add_argument(
        "--grayscale-floor",
        action="store_true",
        default=True,
        help="Use grayscale for floor projections while keeping colorful surfaces.",
    )
    return parser.parse_args()


def _format_time(time_value: float) -> str:
    return f"{time_value:g}"


def _value_limits(field: np.ndarray) -> tuple[float, float]:
    vmin = float(np.nanmin(field))
    vmax = float(np.nanmax(field))
    if np.isclose(vmin, vmax):
        pad = max(abs(vmin), 1.0) * 0.05
        return vmin - pad, vmax + pad
    if vmin >= 0.0:
        return 0.0, vmax + 0.04 * (vmax - vmin)
    pad = 0.04 * (vmax - vmin)
    return vmin - pad, vmax + pad


def _clamp_fraction(value: float, lower: float = 0.12, upper: float = 0.75) -> float:
    return min(max(value, lower), upper)


def _panel_label(label: str, kind: str, time_value: float) -> str:
    kind_label = "Deterministic" if kind == "deterministic" else "Stochastic"
    return f"{kind_label} ($t = {_format_time(time_value)}$)"


def _configure_fonts() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "axes.unicode_minus": False,
        }
    )


def _save_rgb(fig: plt.Figure, output_path: Path, dpi: int) -> None:
    fig.savefig(
        output_path,
        dpi=dpi,
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)
    with Image.open(output_path) as image:
        if image.mode != "RGB":
            image.convert("RGB").save(output_path, dpi=(dpi, dpi))


def _read_endpoint_rows(
    path: Path,
) -> tuple[list[tuple[str, list[dict[str, str]]]], dict[str, str]]:
    metadata: dict[str, str] = {}

    with path.open("r", encoding="utf-8") as handle:
        header: list[str] | None = None
        for raw_line in handle:
            stripped_line = raw_line.strip()
            if stripped_line.startswith("#"):
                key_value = stripped_line[1:].strip().split("=", maxsplit=1)
                if len(key_value) == 2:
                    metadata[key_value[0].strip()] = key_value[1].strip()
                continue
            if stripped_line:
                header = [column.strip() for column in stripped_line.split(",")]
                break

        if header != EXPECTED_SNAPSHOT_HEADER:
            msg = f"Expected CSV header {EXPECTED_SNAPSHOT_HEADER}, got {header}"
            raise ValueError(msg)

        reader = csv.DictReader(handle, fieldnames=header)
        first_key: str | None = None
        first_rows: list[dict[str, str]] = []
        current_key: str | None = None
        current_rows: list[dict[str, str]] = []

        for row in reader:
            time_key = row["time"].strip()
            if current_key is None:
                first_key = time_key
                current_key = time_key
            elif time_key != current_key:
                if current_key == first_key:
                    first_rows = current_rows
                current_key = time_key
                current_rows = []
            current_rows.append(row)

    if current_key is None:
        msg = f"No snapshot rows found in {path}"
        raise ValueError(msg)
    if not first_rows:
        first_rows = current_rows
    return [(first_key or current_key, first_rows), (current_key, current_rows)], metadata


def _rows_to_field(
    rows: list[dict[str, str]],
    path: Path,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    time_value = float(rows[0]["time"])
    nx = max(int(row["x_index"]) for row in rows)
    ny = max(int(row["y_index"]) for row in rows)
    x_coordinates = np.full(nx, np.nan, dtype=float)
    y_coordinates = np.full(ny, np.nan, dtype=float)
    field = np.full((ny, nx), np.nan, dtype=float)

    for row in rows:
        x_index = int(row["x_index"]) - 1
        y_index = int(row["y_index"]) - 1
        x_coordinates[x_index] = float(row["x"])
        y_coordinates[y_index] = float(row["y"])
        field[y_index, x_index] = float(row["value"])

    if np.isnan(x_coordinates).any() or np.isnan(y_coordinates).any():
        msg = f"Incomplete coordinate coverage in {path}"
        raise ValueError(msg)
    if np.isnan(field).any():
        msg = f"Incomplete grid coverage for time={time_value:g} in {path}"
        raise ValueError(msg)

    return time_value, x_coordinates, y_coordinates, field


def _load_endpoint_snapshot_data(
    det_path: Path,
    sto_path: Path,
) -> tuple[
    list[float],
    np.ndarray,
    np.ndarray,
    dict[tuple[str, float], np.ndarray],
    dict[str, str],
]:
    det_endpoint_rows, metadata = _read_endpoint_rows(det_path)
    sto_endpoint_rows, _ = _read_endpoint_rows(sto_path)

    det_snapshots = [_rows_to_field(rows, det_path) for _, rows in det_endpoint_rows]
    sto_snapshots = [_rows_to_field(rows, sto_path) for _, rows in sto_endpoint_rows]
    det_times = [snapshot[0] for snapshot in det_snapshots]
    sto_times = [snapshot[0] for snapshot in sto_snapshots]

    if not np.allclose(det_times, sto_times, rtol=0.0, atol=0.0):
        msg = "Deterministic and stochastic files have mismatched endpoint times."
        raise ValueError(msg)

    x_det = det_snapshots[0][1]
    y_det = det_snapshots[0][2]
    x_sto = sto_snapshots[0][1]
    y_sto = sto_snapshots[0][2]
    if not np.array_equal(x_det, x_sto) or not np.array_equal(y_det, y_sto):
        msg = "Deterministic and stochastic files have mismatched coordinates."
        raise ValueError(msg)

    fields: dict[tuple[str, float], np.ndarray] = {}
    for (time_value, _, _, det_field), (_, _, _, sto_field) in zip(
        det_snapshots,
        sto_snapshots,
    ):
        fields[("deterministic", time_value)] = det_field
        fields[("stochastic", time_value)] = sto_field

    return det_times, x_det, y_det, fields, metadata


def _draw_surface(
    ax: plt.Axes,
    X: np.ndarray,
    Y: np.ndarray,
    field: np.ndarray,
    title: str,
    kind: str,
    args: argparse.Namespace,
    z_limits: tuple[float, float],
    cmap: str,
    panel_id: str | None = None,
    color_limits: tuple[float, float] | None = None,
) -> object:
    z_boost = (
        args.deterministic_z_boost
        if kind == "deterministic"
        else args.stochastic_z_boost
    )
    l_min, l_max = np.nanmin(field), np.nanmax(field)
    if np.isclose(l_min, l_max):
        l_min -= 0.05
        l_max += 0.05
    
    geom_min, geom_max = z_limits
    color_min, color_max = color_limits if color_limits is not None else z_limits
    geom_span = max(geom_max - geom_min, 1.0e-14)
    l_span = max(l_max - l_min, 1.0e-14)
    
    if geom_min < 0:
        z_floor = geom_min - 0.25 * geom_span
    else:
        z_floor = geom_min - 0.15 * geom_span
        
    s = max(1, args.stride)
    
    # Floor: Fine (0.35pt) lines and sparse levels (12 default)
    floor_lw = 0.35 if not args.enhanced_floor else 0.7
    floor_levels = args.contour_levels

    ax.set_proj_type("ortho")
    ax.set_box_aspect((1.0, 1.0, z_boost))
    
    # Surface color may use a row-shared scale, while z geometry is local.
    norm = plt.Normalize(vmin=color_min, vmax=color_max)

    surface = ax.plot_surface(
        X[::s, ::s], Y[::s, ::s], field[::s, ::s],
        cmap=cmap, norm=norm,
        rstride=1, cstride=1, linewidth=0.0,
        antialiased=True, alpha=0.98, shade=True,
    )

    if args.surface_wireframe:
        wire_stride = max(1, args.wireframe_stride)
        ax.plot_wireframe(
            X[::wire_stride, ::wire_stride],
            Y[::wire_stride, ::wire_stride],
            field[::wire_stride, ::wire_stride],
            color="black",
            linewidth=args.wireframe_linewidth,
            alpha=0.32,
        )
    
    # Floor Fill: Grayscale (Greys_r goes from black to white, we want white background)
    # We'll use Greys and ensure transparency or a clean blend.
    l_norm = plt.Normalize(vmin=l_min, vmax=l_max)
    contour_levs = np.linspace(l_min, l_max, floor_levels)
    
    ax.contourf(
        X, Y, field, zdir="z", offset=z_floor,
        levels=contour_levs, cmap="Greys", norm=l_norm, alpha=0.35,
    )
    
    # Floor Borders: Fine BLACK lines (Solid/Dashed)
    pos_levs = [lev for lev in contour_levs if lev > 1e-12]
    neg_levs = [lev for lev in contour_levs if lev < -1e-12]
    
    if pos_levs:
        ax.contour(
            X, Y, field, zdir="z", offset=z_floor,
            levels=pos_levs, colors="black", linewidths=floor_lw,
            alpha=0.8, linestyles="solid"
        )
    if neg_levs:
        ax.contour(
            X, Y, field, zdir="z", offset=z_floor,
            levels=neg_levs, colors="black", linewidths=floor_lw,
            alpha=0.8, linestyles="dashed"
        )

    # Surface contours: Draw fine lines on the mesh to show structure
    ax.contour(
        X, Y, field,
        levels=contour_levs, colors="black", linewidths=0.25,
        alpha=0.4
    )

    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlim(float(np.nanmin(X)), float(np.nanmax(X)))
    ax.set_ylim(float(np.nanmin(Y)), float(np.nanmax(Y)))
    ax.set_zlim(z_floor, l_max + 0.05 * l_span)
    
    ax.set_xlabel("$x$", fontsize=9, labelpad=4)
    ax.set_ylabel("$y$", fontsize=9, labelpad=4)
    ax.set_zlabel("")
    ax.set_title(title, fontsize=10, pad=3)
    
    from matplotlib import ticker
    ax.tick_params(axis="both", labelsize=7.5, pad=2)
    ax.zaxis.set_tick_params(labelsize=7.5, pad=2)
    if l_span < 0.01:
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter("%.1e"))
    elif l_span < 1.0:
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    else:
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        
    ax.locator_params(axis="x", nbins=4)
    ax.locator_params(axis="y", nbins=4)
    ax.locator_params(axis="z", nbins=4)

    if panel_id:
        ax.text2D(-0.10, 1.04, panel_id, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="right")

    return surface


def _zoom_window(field: np.ndarray, fraction: float) -> tuple[slice, slice]:
    n_rows, n_cols = field.shape
    row_center, col_center = np.unravel_index(np.nanargmax(np.abs(field)), field.shape)
    row_span = max(3, int(round(_clamp_fraction(fraction) * n_rows)))
    col_span = max(3, int(round(_clamp_fraction(fraction) * n_cols)))
    row_start = min(max(row_center - row_span // 2, 0), max(n_rows - row_span, 0))
    col_start = min(max(col_center - col_span // 2, 0), max(n_cols - col_span, 0))
    return slice(row_start, row_start + row_span), slice(col_start, col_start + col_span)


def _add_panel_inset(
    fig: plt.Figure,
    ax: plt.Axes,
    X: np.ndarray,
    Y: np.ndarray,
    field: np.ndarray,
    cmap: str,
    args: argparse.Namespace,
    is_combined: bool = False,
) -> None:
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    row_slice, col_slice = _zoom_window(field, args.inset_zoom_fraction)
    X_zoom = X[row_slice, col_slice]
    Y_zoom = Y[row_slice, col_slice]
    field_zoom = field[row_slice, col_slice]
    zoom_min, zoom_max = _value_limits(field_zoom)
    width = "28%" if is_combined else "30%"
    height = "28%" if is_combined else "30%"
    ax_ins = inset_axes(ax, width=width, height=height, loc="upper right", borderpad=1.15)
    ax_ins.imshow(
        field_zoom,
        cmap=cmap,
        origin="lower",
        interpolation="bilinear",
        aspect="auto",
        vmin=zoom_min,
        vmax=zoom_max,
        extent=[
            float(np.nanmin(X_zoom)),
            float(np.nanmax(X_zoom)),
            float(np.nanmin(Y_zoom)),
            float(np.nanmax(Y_zoom)),
        ],
    )
    if not np.isclose(np.nanmin(field_zoom), np.nanmax(field_zoom)):
        inset_levels = np.linspace(np.nanmin(field_zoom), np.nanmax(field_zoom), 7)
        ax_ins.contour(
            X_zoom,
            Y_zoom,
            field_zoom,
            levels=inset_levels,
            colors="black",
            linewidths=0.35,
            alpha=0.65,
        )
    ax_ins.set_xticks([])
    ax_ins.set_yticks([])
    for spine in ax_ins.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.8)
    ax_ins.set_title("local zoom", fontsize=7, pad=2, fontweight="bold")


def _draw_panel(
    output_path: Path,
    X: np.ndarray,
    Y: np.ndarray,
    field: np.ndarray,
    title: str,
    kind: str,
    args: argparse.Namespace,
    z_limits: tuple[float, float],
    cmap: str,
) -> None:
    fig = plt.figure(figsize=(args.width_inches, args.height_inches), facecolor="white", constrained_layout=False)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    surface = _draw_surface(ax, X, Y, field, title, kind, args, z_limits=z_limits, cmap=cmap)
    if args.insets:
        _add_panel_inset(fig, ax, X, Y, field, cmap, args, is_combined=False)
    colorbar = fig.colorbar(surface, ax=ax, fraction=0.036, pad=0.04, shrink=0.75)
    colorbar.ax.tick_params(labelsize=7.5, width=0.6, pad=2)
    colorbar.set_label("Solution $u$", fontsize=9, labelpad=3)
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.12, top=0.92)
    _save_rgb(fig, output_path, args.dpi)


def _draw_combined(
    output_path: Path,
    X: np.ndarray,
    Y: np.ndarray,
    panel_fields: list[tuple[str, str, float, np.ndarray]],
    variant_name: str,
    args: argparse.Namespace,
    z_limits_by_kind: dict[str, tuple[float, float]],
    cmaps_by_kind: dict[str, str],
) -> None:
    fig = plt.figure(figsize=(7.2, 6.8), facecolor="white", constrained_layout=False)
    surfaces = {}
    for panel_index, (label, kind, time_value, field) in enumerate(panel_fields):
        ax = fig.add_subplot(2, 2, panel_index + 1, projection="3d")
        panel_id = PANEL_IDS[panel_index]
        cmap = cmaps_by_kind[kind]
        surface = _draw_surface(
            ax,
            X,
            Y,
            field,
            _panel_label(label, kind, time_value),
            kind,
            args,
            z_limits=_value_limits(field),
            cmap=cmap,
            panel_id=panel_id,
            color_limits=z_limits_by_kind[kind],
        )
        surfaces[kind] = surface
        if args.insets:
            _add_panel_inset(fig, ax, X, Y, field, cmap, args, is_combined=True)

    cax_det = fig.add_axes([0.90, 0.58, 0.02, 0.32])
    cb_det = fig.colorbar(surfaces["deterministic"], cax=cax_det)
    cb_det.ax.tick_params(labelsize=7, width=0.6, pad=2)
    cb_det.set_label("Det. row scale", fontsize=8, labelpad=3)
    
    cax_sto = fig.add_axes([0.90, 0.12, 0.02, 0.32])
    cb_sto = fig.colorbar(surfaces["stochastic"], cax=cax_sto)
    cb_sto.ax.tick_params(labelsize=7, width=0.6, pad=2)
    cb_sto.set_label("Sto. row scale", fontsize=8, labelpad=3)

    fig.subplots_adjust(left=0.08, right=0.87, bottom=0.08, top=0.94, wspace=0.14, hspace=0.22)
    _save_rgb(fig, output_path, args.dpi)


def _variant_names(argument: str) -> list[str]:
    if argument == "both": return ["signed", "nonnegative"]
    return [argument]

def _layout_names(argument: str) -> list[str]:
    if argument == "both": return ["separate", "combined"]
    return [argument]

def _variant_field(field: np.ndarray, variant_name: str) -> np.ndarray:
    if variant_name == "nonnegative": return np.maximum(field, 0.0)
    return field

def _get_z_limits_by_kind(panel_fields: list[tuple[str, str, float, np.ndarray]]) -> dict[str, tuple[float, float]]:
    limits = {}
    for kind in ["deterministic", "stochastic"]:
        vals = np.concatenate([f.ravel() for _, k, _, f in panel_fields if k == kind])
        limits[kind] = _value_limits(vals)
    return limits

def main() -> None:
    args = parse_args()
    if args.deterministic_inset:
        args.insets = True
    args.inset_zoom_fraction = _clamp_fraction(args.inset_zoom_fraction)
    args.wireframe_stride = max(1, args.wireframe_stride)
    args.wireframe_linewidth = max(0.0, args.wireframe_linewidth)
    det_path, sto_path = resolve_input_paths(args.input)
    prefix = _extract_timestamp_prefix(det_path)
    output_dir = Path(args.output_dir) if args.output_dir else DEFAULT_PLOT_DIR / f"{prefix}_manuscript_panels"
    output_dir.mkdir(parents=True, exist_ok=True)
    times, x_coords, y_coords, fields, _ = _load_endpoint_snapshot_data(det_path, sto_path)
    X, Y = np.meshgrid(x_coords, y_coords)
    _configure_fonts()
    requested_layouts = set(_layout_names(args.layout))
    written_paths: list[Path] = []
    
    if args.cmap: cmaps_by_kind = {"deterministic": args.cmap, "stochastic": args.cmap}
    else: cmaps_by_kind = {"deterministic": args.cmap_deterministic, "stochastic": args.cmap_stochastic}
        
    for variant_name in _variant_names(args.variant):
        variant_dir = output_dir / variant_name
        panel_fields: list[tuple[str, str, float, np.ndarray]] = []
        for label, kind, time_index in PANEL_SPECS:
            time_value = times[time_index]
            field = _variant_field(np.array(fields[(kind, time_value)], copy=True), variant_name)
            panel_fields.append((label, kind, time_value, field))

        z_limits_by_kind = _get_z_limits_by_kind(panel_fields)

        if "separate" in requested_layouts:
            separate_dir = variant_dir / "separate"
            separate_dir.mkdir(parents=True, exist_ok=True)
            for label, kind, time_value, field in panel_fields:
                output_path = separate_dir / f"{prefix}_{variant_name}_{label}_{kind}.{args.extension}"
                local_limits = _value_limits(field)
                _draw_panel(output_path, X, Y, field, _panel_label(label, kind, time_value), kind, args, z_limits=local_limits, cmap=cmaps_by_kind[kind])
                written_paths.append(output_path)

        if "combined" in requested_layouts:
            variant_dir.mkdir(parents=True, exist_ok=True)
            output_path = variant_dir / f"{prefix}_{variant_name}_combined_2x2.{args.extension}"
            _draw_combined(output_path, X, Y, panel_fields, variant_name, args, z_limits_by_kind=z_limits_by_kind, cmaps_by_kind=cmaps_by_kind)
            written_paths.append(output_path)

    print("Wrote manuscript panels:")
    for path in written_paths: print(f"  {path}")

if __name__ == "__main__":
    main()
