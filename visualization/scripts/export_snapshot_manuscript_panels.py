#!/usr/bin/env python3
"""Export manuscript-ready 3-D snapshot panels."""

from __future__ import annotations

import argparse
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
    load_snapshot_data,
    resolve_input_paths,
)


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
        help="Add a local-scale inset to deterministic panels to show structure.",
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
    
    zmin, zmax = z_limits
    z_span = max(zmax - zmin, 1.0e-14)
    l_span = max(l_max - l_min, 1.0e-14)
    
    if zmin < 0:
        z_floor = zmin - 0.25 * z_span
    else:
        z_floor = zmin - 0.15 * z_span
        
    s = max(1, args.stride)
    
    # Floor: Fine (0.35pt) lines and sparse levels (12 default)
    floor_lw = 0.35 if not args.enhanced_floor else 0.7
    floor_levels = args.contour_levels

    ax.set_proj_type("ortho")
    ax.set_box_aspect((1.0, 1.0, z_boost))
    
    # Surface: Colorful (turbo/Spectral_r)
    norm = plt.Normalize(vmin=zmin, vmax=zmax)

    surface = ax.plot_surface(
        X[::s, ::s], Y[::s, ::s], field[::s, ::s],
        cmap=cmap, norm=norm,
        rstride=1, cstride=1, linewidth=0.0,
        antialiased=True, alpha=0.98, shade=True,
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


def _add_deterministic_inset(
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
    ax_ins = inset_axes(ax, width="25%", height="25%", loc="upper right", borderpad=1.2)
    ax_ins.imshow(field, cmap=cmap, origin="lower", interpolation="bilinear", aspect="auto")
    ax_ins.set_xticks([])
    ax_ins.set_yticks([])
    for spine in ax_ins.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.8)
    ax_ins.set_title("zoom (local)", fontsize=7, pad=2, fontweight="bold")


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
    if args.deterministic_inset and kind == "deterministic":
        _add_deterministic_inset(fig, ax, X, Y, field, cmap, args, is_combined=False)
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
    fig = plt.figure(figsize=(7.8, 7.2), facecolor="white", constrained_layout=False)
    surfaces = {}
    for panel_index, (label, kind, time_value, field) in enumerate(panel_fields):
        ax = fig.add_subplot(2, 2, panel_index + 1, projection="3d")
        panel_id = PANEL_IDS[panel_index]
        cmap = cmaps_by_kind[kind]
        surface = _draw_surface(ax, X, Y, field, _panel_label(label, kind, time_value), kind, args, z_limits=z_limits_by_kind[kind], cmap=cmap, panel_id=panel_id)
        surfaces[kind] = surface
        if args.deterministic_inset and kind == "deterministic":
            _add_deterministic_inset(fig, ax, X, Y, field, cmap, args, is_combined=True)

    cax_det = fig.add_axes([0.91, 0.58, 0.02, 0.32])
    cb_det = fig.colorbar(surfaces["deterministic"], cax=cax_det)
    cb_det.ax.tick_params(labelsize=7, width=0.6, pad=2)
    cb_det.set_label("Det. row scale", fontsize=8, labelpad=3)
    
    cax_sto = fig.add_axes([0.91, 0.12, 0.02, 0.32])
    cb_sto = fig.colorbar(surfaces["stochastic"], cax=cax_sto)
    cb_sto.ax.tick_params(labelsize=7, width=0.6, pad=2)
    cb_sto.set_label("Sto. row scale", fontsize=8, labelpad=3)

    fig.subplots_adjust(left=0.10, right=0.88, bottom=0.10, top=0.94, wspace=0.16, hspace=0.25)
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
    det_path, sto_path = resolve_input_paths(args.input)
    prefix = _extract_timestamp_prefix(det_path)
    output_dir = Path(args.output_dir) if args.output_dir else DEFAULT_PLOT_DIR / f"{prefix}_manuscript_panels"
    output_dir.mkdir(parents=True, exist_ok=True)
    times, x_coords, y_coords, fields, _ = load_snapshot_data(det_path, sto_path)
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
