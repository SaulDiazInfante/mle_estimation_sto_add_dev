#!/usr/bin/env python3
"""3D surface comparison: deterministic vs stochastic nonnegative solution.

Produces two side-by-side 3D surface panels:
  1. Deterministic  – smooth reference surface
  2. Stochastic     – same field perturbed by noise (visually rugose)

Usage
-----
    python plot_3d_rugosity_comparison.py            # latest CSV, last snapshot
    python plot_3d_rugosity_comparison.py --time-index 50
    python plot_3d_rugosity_comparison.py --elev 35 --azim 210 --stride 2
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 – registers 3-D projection

# Reuse loaders/helpers from the 2-D companion script
_SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(_SCRIPTS_DIR))
from plot_solution_snapshot_comparison import (  # noqa: E402
    DEFAULT_PLOT_DIR,
    _extract_timestamp_prefix,
    load_snapshot_data,
    resolve_input_paths,
)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="3-D comparison: deterministic (smooth) vs stochastic (rough)."
    )
    p.add_argument(
        "--input", default=None,
        help="Snapshot CSV; defaults to the latest timestamped file.",
    )
    p.add_argument(
        "--output", default=None,
        help="Output PNG path; defaults to <stem>_3d_rugosity.png beside the input.",
    )
    p.add_argument(
        "--time-index", type=int, default=-1,
        help="Index of the time snapshot to visualise (-1 = last, most evolved).",
    )
    p.add_argument(
        "--elev", type=float, default=28.0,
        help="Camera elevation angle in degrees.",
    )
    p.add_argument(
        "--azim", type=float, default=225.0,
        help="Camera azimuth angle in degrees.",
    )
    p.add_argument(
        "--stride", type=int, default=1,
        help="Mesh stride (1 = full resolution; increase to speed up rendering).",
    )
    p.add_argument(
        "--alpha", type=float, default=0.92,
        help="Surface opacity [0..1].",
    )
    p.add_argument(
        "--cmap", default="magma",
        help="Colormap for the nonnegative solution surfaces.",
    )
    p.add_argument(
        "--contour-levels", type=int, default=18,
        help="Number of filled contour levels projected onto the xy plane.",
    )
    p.add_argument(
        "--dpi", type=int, default=200,
        help="Output image resolution.",
    )
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def _fmt_t(t: float) -> str:
    return f"{t:g}"


def _build_surface(
    ax: plt.Axes,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    cmap_name: str,
    vmin: float,
    vmax: float,
    stride: int,
    alpha: float,
    elev: float,
    azim: float,
    title: str,
    zlabel: str,
    contour_levels: int,
) -> None:
    s = stride
    z_span = max(abs(vmax - vmin), 1e-14)
    z_floor = 0.0
    ax.plot_surface(
        X[::s, ::s],
        Y[::s, ::s],
        Z[::s, ::s],
        cmap=cmap_name,
        vmin=vmin,
        vmax=vmax,
        rstride=1,
        cstride=1,
        antialiased=False,
        linewidth=0.0,
        alpha=alpha,
    )
    ax.contourf(
        X,
        Y,
        Z,
        zdir="z",
        offset=z_floor,
        levels=contour_levels,
        cmap=cmap_name,
        vmin=vmin,
        vmax=vmax,
        alpha=0.80,
    )
    ax.contour(
        X,
        Y,
        Z,
        zdir="z",
        offset=z_floor,
        levels=max(6, contour_levels // 2),
        colors="white",
        linewidths=0.45,
        alpha=0.60,
    )
    ax.plot_wireframe(
        X[::s, ::s],
        Y[::s, ::s],
        Z[::s, ::s],
        rstride=max(1, 4 * s),
        cstride=max(1, 4 * s),
        color="k",
        linewidth=0.16,
        alpha=0.26,
    )
    z_pad = 0.06 * z_span
    ax.set_zlim(z_floor, vmax + z_pad)
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(title, fontsize=11, pad=6)
    ax.set_xlabel(r"$x$", labelpad=4)
    ax.set_ylabel(r"$y$", labelpad=4)
    ax.set_zlabel(zlabel, labelpad=4)
    ax.tick_params(axis="both", labelsize=7, pad=2)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    if args.stride < 1:
        raise ValueError("--stride must be at least 1.")
    if args.contour_levels < 2:
        raise ValueError("--contour-levels must be at least 2.")

    det_path, sto_path = resolve_input_paths(args.input)

    times, x_coords, y_coords, fields, _ = load_snapshot_data(det_path, sto_path)

    t_idx = args.time_index % len(times)
    t = times[t_idx]

    det = np.maximum(fields[("deterministic", t)], 0.0)
    sto = np.maximum(fields[("stochastic", t)], 0.0)

    X, Y = np.meshgrid(x_coords, y_coords)

    # Shared colour / z limits for the det / sto panels after physical masking.
    vmin = 0.0
    vmax = max(det.max(), sto.max())
    if vmax <= vmin:
        vmax = vmin + 1.0e-14

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(13.2, 6.4))

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    _build_surface(
        ax1, X, Y, det,
        cmap_name=args.cmap,
        vmin=vmin, vmax=vmax,
        stride=args.stride, alpha=args.alpha,
        elev=args.elev, azim=args.azim,
        title="Deterministic  (smooth)",
        zlabel=r"$u_{\mathrm{det}}$",
        contour_levels=args.contour_levels,
    )

    _build_surface(
        ax2, X, Y, sto,
        cmap_name=args.cmap,
        vmin=vmin, vmax=vmax,
        stride=args.stride, alpha=args.alpha,
        elev=args.elev, azim=args.azim,
        title="Stochastic  (rugose)",
        zlabel=r"$u_{\mathrm{sto}}$",
        contour_levels=args.contour_levels,
    )

    # ── Colorbars ─────────────────────────────────────────────────────────────
    # Shared bar for det + sto
    sm_uv = plt.cm.ScalarMappable(
        cmap=args.cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax)
    )
    sm_uv.set_array([])
    cb1 = fig.colorbar(sm_uv, ax=[ax1, ax2], fraction=0.022, pad=0.06, shrink=0.55)
    cb1.set_label(r"$u$ (negative values masked to 0)", fontsize=10)
    cb1.ax.tick_params(labelsize=8)

    # ── Title & layout ────────────────────────────────────────────────────────
    n_times = len(times)
    fig.suptitle(
        rf"3-D nonnegative solution comparison   $t = {_fmt_t(t)}$"
        rf"   (snapshot {t_idx % n_times + 1}/{n_times})",
        fontsize=13,
        y=1.00,
    )
    fig.subplots_adjust(left=0.03, right=0.88, bottom=0.05, top=0.94, wspace=0.04)

    # ── Save ──────────────────────────────────────────────────────────────────
    if args.output is not None:
        output_path = Path(args.output)
    else:
        prefix = _extract_timestamp_prefix(det_path)
        output_path = DEFAULT_PLOT_DIR / f"{prefix}_3d_rugosity.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote 3-D rugosity plot to {output_path}")


if __name__ == "__main__":
    main()
