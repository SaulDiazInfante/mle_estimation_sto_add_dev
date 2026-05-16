#!/usr/bin/env python3
"""3D surface comparison: smooth deterministic vs rugose stochastic solution.

Produces three side-by-side 3D surface panels:
  1. Deterministic  – smooth reference surface
  2. Stochastic     – same field perturbed by noise (visually rugose)
  3. Noise effect   – pointwise difference  u_sto - u_det  (signed, diverging palette)

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
    load_snapshot_data,
    resolve_input_path,
)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="3-D rugosity comparison: deterministic (smooth) vs stochastic (rough)."
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
        "--dpi", type=int, default=200,
        help="Output image resolution.",
    )
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def _fmt_t(t: float) -> str:
    return f"{t:g}"


def _safe_norm(arr: np.ndarray, lo: float, hi: float) -> np.ndarray:
    span = hi - lo
    if span < 1e-14:
        return np.full_like(arr, 0.5)
    return np.clip((arr - lo) / span, 0.0, 1.0)


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
) -> None:
    s = stride
    surf = ax.plot_surface(
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
    # Subtle wireframe outline for the rugosity panels
    if cmap_name != "RdBu_r":
        ax.plot_wireframe(
            X[::s, ::s],
            Y[::s, ::s],
            Z[::s, ::s],
            rstride=max(1, 4 * s),
            cstride=max(1, 4 * s),
            color="k",
            linewidth=0.18,
            alpha=0.30,
        )
    z_pad = 0.06 * max(abs(vmax - vmin), 1e-14)
    ax.set_zlim(vmin - z_pad, vmax + z_pad)
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(title, fontsize=11, pad=6)
    ax.set_xlabel(r"$x$", labelpad=4)
    ax.set_ylabel(r"$y$", labelpad=4)
    ax.set_zlabel(zlabel, labelpad=4)
    ax.tick_params(axis="both", labelsize=7, pad=2)
    return surf


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)

    times, x_coords, y_coords, fields, _ = load_snapshot_data(input_path)

    t_idx = args.time_index % len(times)
    t = times[t_idx]

    det = fields[("deterministic", t)]
    sto = fields[("stochastic", t)]
    diff = sto - det

    X, Y = np.meshgrid(x_coords, y_coords)

    # Shared colour / z limits for the det / sto panels
    vmin = min(det.min(), sto.min())
    vmax = max(det.max(), sto.max())

    # Symmetric limits for the difference panel
    d_lim = max(abs(diff.min()), abs(diff.max()), 1e-14)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(19, 6.4))

    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")

    surf_det = _build_surface(
        ax1, X, Y, det,
        cmap_name="viridis",
        vmin=vmin, vmax=vmax,
        stride=args.stride, alpha=args.alpha,
        elev=args.elev, azim=args.azim,
        title="Deterministic  (smooth)",
        zlabel=r"$u_{\mathrm{det}}$",
    )

    surf_sto = _build_surface(
        ax2, X, Y, sto,
        cmap_name="viridis",
        vmin=vmin, vmax=vmax,
        stride=args.stride, alpha=args.alpha,
        elev=args.elev, azim=args.azim,
        title="Stochastic  (rugose)",
        zlabel=r"$u_{\mathrm{sto}}$",
    )

    surf_diff = _build_surface(
        ax3, X, Y, diff,
        cmap_name="RdBu_r",
        vmin=-d_lim, vmax=d_lim,
        stride=args.stride, alpha=args.alpha,
        elev=args.elev, azim=args.azim,
        title=r"Noise effect  ($u_{\mathrm{sto}} - u_{\mathrm{det}}$)",
        zlabel=r"$\Delta u$",
    )

    # ── Colorbars ─────────────────────────────────────────────────────────────
    # Shared bar for det + sto
    sm_uv = plt.cm.ScalarMappable(
        cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax)
    )
    sm_uv.set_array([])
    cb1 = fig.colorbar(sm_uv, ax=[ax1, ax2], fraction=0.022, pad=0.06, shrink=0.55)
    cb1.set_label(r"$u$", fontsize=10)
    cb1.ax.tick_params(labelsize=8)

    sm_rb = plt.cm.ScalarMappable(
        cmap="RdBu_r", norm=plt.Normalize(vmin=-d_lim, vmax=d_lim)
    )
    sm_rb.set_array([])
    cb3 = fig.colorbar(sm_rb, ax=ax3, fraction=0.030, pad=0.06, shrink=0.55)
    cb3.set_label(r"$\Delta u$", fontsize=10)
    cb3.ax.tick_params(labelsize=8)

    # ── Title & layout ────────────────────────────────────────────────────────
    n_times = len(times)
    fig.suptitle(
        rf"3-D rugosity comparison   $t = {_fmt_t(t)}$"
        rf"   (snapshot {t_idx % n_times + 1}/{n_times})",
        fontsize=13,
        y=1.00,
    )
    fig.subplots_adjust(left=0.02, right=0.88, bottom=0.04, top=0.94, wspace=0.05)

    # ── Save ──────────────────────────────────────────────────────────────────
    if args.output is not None:
        output_path = Path(args.output)
    else:
        output_path = DEFAULT_PLOT_DIR / f"{input_path.stem}_3d_rugosity.png"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote 3-D rugosity plot to {output_path}")


if __name__ == "__main__":
    main()
