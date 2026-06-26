#!/usr/bin/env python3
"""Plot the time evolution of sigma_hat, beta_hat, and theta_hat."""

from __future__ import annotations

import argparse
import csv
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize, SymLogNorm, PowerNorm
from matplotlib.patches import Patch

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PLOT_DIR = Path("visualization/plots")
ESTIMATOR_SUFFIX = "_estimator_trajectory.csv"
PARAMETER_PANELS = (
    {
        "estimate_key": "sigma_hat",
        "truth_key": "sigma_true",
        "symbol": r"\sigma",
        "panel_label": "(A)",
    },
    {
        "estimate_key": "beta_hat",
        "truth_key": "beta_true",
        "symbol": r"\beta",
        "panel_label": "(B)",
    },
    {
        "estimate_key": "theta_hat",
        "truth_key": "theta_true",
        "symbol": r"\theta",
        "panel_label": "(C)",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot estimator trajectories against their true values."
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Input trajectory file; defaults to the latest timestamped estimator CSV.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output image path; defaults to a matching timestamped PNG in visualization/plots.",
    )
    parser.add_argument(
        "--x-axis",
        choices=("n_obs", "time"),
        default="n_obs",
        help="Horizontal axis for the trajectory plot.",
    )
    parser.add_argument(
        "--mode",
        choices=("value", "normalized", "relative-error"),
        default="value",
        help="Quantity to plot for each parameter trajectory.",
    )
    parser.add_argument(
        "--y-scale",
        choices=("linear", "log", "auto"),
        default="log",
        help="Vertical scaling. Use auto to restore sign-aware automatic scaling.",
    )
    parser.add_argument(
        "--tolerance-band",
        type=float,
        default=0.003,
        help="Fractional tolerance band around the reference value; set to 0 to disable.",
    )
    parser.add_argument(
        "--min-observation",
        type=int,
        default=None,
        help="Minimum n_obs row to include in the plot; defaults to the first n_obs in the input.",
    )
    parser.add_argument(
        "--max-observation",
        type=int,
        default=40000,
        help="Maximum n_obs row to include in the plot.",
    )
    parser.add_argument(
        "--no-tolerance-background",
        dest="tolerance_background",
        action="store_false",
        help="Disable the tolerance-distance background colormap.",
    )
    parser.add_argument(
        "--tolerance-background-cmap",
        default="Spectral",
        help="Matplotlib colormap for the tolerance-distance background.",
    )
    parser.add_argument(
        "--tolerance-background-max-multiple",
        type=float,
        default=0.2,
        help="Tolerance multiple mapped to the farthest background color.",
    )
    parser.add_argument(
        "--tolerance-background-scale",
        choices=("linear", "log", "power"),
        default="power",
        help="Color normalization for the tolerance-distance background.",
    )
    parser.add_argument(
        "--tolerance-background-alpha",
        type=float,
        default=0.34,
        help="Opacity for the tolerance-distance background.",
    )
    parser.add_argument(
        "--tolerance-background-normalization",
        choices=("panel", "absolute"),
        default="panel",
        help="Normalize colors per panel or against the absolute max multiple.",
    )
    parser.set_defaults(tolerance_background=True)
    return parser.parse_args()


def resolve_input_path(argument: str | None) -> Path:
    if argument is not None:
        return Path(argument)

    candidates = sorted(DEFAULT_INPUT_DIR.glob(f"*{ESTIMATOR_SUFFIX}"))
    if not candidates:
        msg = (
            "No generated estimator data was found in "
            f"{DEFAULT_INPUT_DIR}. First run the application with 'make run'."
        )
        raise FileNotFoundError(msg)

    print(f"Using latest estimator history from {candidates[-1]}")
    return candidates[-1]


def resolve_output_path(argument: str | None, input_path: Path) -> Path:
    if argument is not None:
        return Path(argument)

    if input_path.name.endswith(".csv"):
        output_name = f"{input_path.stem}.png"
    else:
        output_name = f"{input_path.name}.png"

    return DEFAULT_PLOT_DIR / output_name


def load_trajectory(path: Path) -> dict[str, list[float]]:
    data = {
        "n_obs": [],
        "time": [],
        "sigma_hat": [],
        "sigma_true": [],
        "beta_hat": [],
        "beta_true": [],
        "theta_hat": [],
        "theta_true": [],
    }

    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = list(data.keys())
        if reader.fieldnames != expected:
            msg = f"Expected CSV header {expected}, got {reader.fieldnames}"
            raise ValueError(msg)

        for row in reader:
            data["n_obs"].append(int(row["n_obs"]))
            data["time"].append(float(row["time"]))
            data["sigma_hat"].append(float(row["sigma_hat"]))
            data["sigma_true"].append(float(row["sigma_true"]))
            data["beta_hat"].append(float(row["beta_hat"]))
            data["beta_true"].append(float(row["beta_true"]))
            data["theta_hat"].append(float(row["theta_hat"]))
            data["theta_true"].append(float(row["theta_true"]))

    if not data["time"]:
        msg = f"No trajectory rows found in {path}"
        raise ValueError(msg)

    return data


def filter_observation_window(
    trajectory: dict[str, list[float]],
    min_observation: int,
    max_observation: int,
) -> dict[str, list[float]]:
    if min_observation > max_observation:
        msg = (
            "Minimum observation cannot exceed maximum observation: "
            f"{min_observation} > {max_observation}"
        )
        raise ValueError(msg)

    selected_indices = [
        index
        for index, n_obs in enumerate(trajectory["n_obs"])
        if min_observation <= n_obs <= max_observation
    ]
    if not selected_indices:
        msg = (
            "No trajectory rows fall inside the requested observation window "
            f"[{min_observation}, {max_observation}]"
        )
        raise ValueError(msg)

    return {
        key: [values[index] for index in selected_indices]
        for key, values in trajectory.items()
    }


def configure_y_scale(
    ax: plt.Axes,
    estimate: list[float],
    truth: list[float],
    title: str,
    y_scale: str,
) -> None:
    if y_scale == "linear":
        ax.set_yscale("linear")
        return

    if y_scale == "log":
        ax.set_yscale("log")
        return

    values = estimate + truth
    nonzero_values = [abs(value) for value in values if value != 0.0]

    if not nonzero_values:
        ax.set_yscale("linear")
        return

    min_abs_nonzero = min(nonzero_values)

    if all(value > 0.0 for value in values):
        ax.set_yscale("log")
        return

    if all(value < 0.0 for value in values):
        ax.set_yscale("log")
        ax.invert_yaxis()
        return

    # Fall back to a symmetric log scale if the estimator changes sign.
    ax.set_yscale("symlog", linthresh=max(min_abs_nonzero, 1.0e-12))
    print(
        f"Using symlog on the {title} panel because the series is not strictly positive."
    )


def choose_x_axis(
    trajectory: dict[str, list[float]],
    x_axis: str,
) -> tuple[Sequence[float] | Sequence[int], str]:
    if x_axis == "time":
        return trajectory["time"], r"Time $t$"

    return trajectory["n_obs"], "Number of observations"


def transform_series(
    estimate: list[float],
    truth: list[float],
    mode: str,
) -> tuple[list[float], list[float]]:
    if mode == "value":
        return estimate, truth

    transformed_estimate: list[float] = []
    transformed_truth: list[float] = []

    for estimate_value, truth_value in zip(estimate, truth):
        if truth_value == 0.0:
            signed_scale = 1.0
        else:
            signed_scale = truth_value
        absolute_scale = abs(signed_scale)

        if mode == "normalized":
            transformed_estimate.append(estimate_value / signed_scale)
            transformed_truth.append(truth_value / signed_scale)
        else:
            transformed_estimate.append(abs(estimate_value - truth_value) / absolute_scale)
            transformed_truth.append(0.0)

    return transformed_estimate, transformed_truth


def y_axis_label(mode: str, symbol: str) -> str:
    if mode == "normalized":
        return "Estimated / true value"
    if mode == "relative-error":
        return "Relative error"
    return rf"$Estimated\ \widehat{{{symbol}}}$ : "


def reference_label(mode: str) -> str:
    if mode == "relative-error":
        return "Zero error"
    if mode == "normalized":
        return "Reference"
    return "True value"


def format_tolerance_label(tolerance_band: float) -> str:
    percent_value = tolerance_band * 100.0
    if percent_value >= 1.0:
        return f"{percent_value:.0f}% tolerance band"
    return f"{percent_value:.2f}% tolerance band"


def build_tolerance_norm(max_rel_error: float, scale: str, vmin: float = 0.0) -> Normalize:
    if scale == "log":
        return SymLogNorm(
            linthresh=max_rel_error * 0.01,
            vmin=vmin,
            vmax=max_rel_error,
            base=10.0,
            clip=True,
        )
    if scale == "power":
        return PowerNorm(
            gamma=0.5,
            vmin=vmin,
            vmax=max_rel_error,
            clip=True,
        )

    return Normalize(vmin=vmin, vmax=max_rel_error, clip=True)


def tolerance_colorbar_ticks(max_rel_error: float, scale: str, vmin: float = 0.0) -> list[float]:
    if scale == "log":
        return [vmin, max_rel_error * 0.1, max_rel_error]

    return [vmin, 0.5 * max_rel_error, max_rel_error]


def tolerance_distance_multiple(
    value: float,
    reference_value: float,
    mode: str,
    tolerance_band: float,
) -> float:
    if mode == "relative-error":
        fractional_distance = abs(value)
    else:
        reference_scale = max(abs(reference_value), 1.0e-12)
        fractional_distance = abs(value - reference_value) / reference_scale

    return fractional_distance / tolerance_band


def quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        return 0.0

    sorted_values = sorted(values)
    index = probability * (len(sorted_values) - 1)
    lower_index = int(index)
    upper_index = min(lower_index + 1, len(sorted_values) - 1)
    fraction = index - lower_index
    return (
        sorted_values[lower_index] * (1.0 - fraction)
        + sorted_values[upper_index] * fraction
    )


def panel_tolerance_norm(multiples: Sequence[float], scale: str) -> Normalize:
    low = quantile(multiples, 0.05)
    high = quantile(multiples, 0.95)
    if high <= low:
        high = max(low + 1.0e-12, low * 1.001)

    if scale == "log":
        return SymLogNorm(
            linthresh=max((high - low) * 0.01, 1.0e-12),
            vmin=low,
            vmax=high,
            base=10.0,
            clip=True,
        )

    return Normalize(vmin=low, vmax=high, clip=True)


def x_value_edges(x_values: Sequence[float] | Sequence[int]) -> list[float]:
    values = [float(value) for value in x_values]
    if len(values) == 1:
        return [values[0] - 0.5, values[0] + 0.5]

    edges = [values[0] - 0.5 * (values[1] - values[0])]
    edges.extend(
        0.5 * (left + right)
        for left, right in zip(values[:-1], values[1:])
    )
    edges.append(values[-1] + 0.5 * (values[-1] - values[-2]))
    return edges


def plot_tolerance_background(
    ax: plt.Axes,
    x_values: Sequence[float] | Sequence[int],
    estimate: list[float],
    reference: list[float],
    mode: str,
    tolerance_band: float,
    cmap: Colormap,
    norm: Normalize,
    normalization: str,
    scale: str,
    alpha: float,
) -> None:
    if tolerance_band <= 0.0:
        return

    lower, upper = ax.get_ylim()
    if lower == upper:
        return

    # Create a grid of values for the background
    # The color depends only on the y-value (the distance from the truth)
    import numpy as np
    y_res = 200
    y_vals = np.linspace(lower, upper, y_res)
    
    # Calculate the relative error for each y_val
    # truth is constant for a given panel
    truth = reference[0]
    
    # Background color as a function of y only
    rel_errors = []
    for y in y_vals:
        # relative error = |y - truth| / truth * 100
        rel_err = (abs(y - truth) / max(abs(truth), 1.0e-12)) * 100.0
        rel_errors.append(rel_err)
    
    # To use pcolormesh, we need coordinates of the corners (edges).
    # For a single strip across X:
    xmin, xmax = ax.get_xlim()
    X = np.array([xmin, xmax])
    
    # Y contains the edges of the horizontal strips.
    # Y has length y_res.
    Y = y_vals
    
    # C contains the colors for each strip.
    # C has shape (y_res - 1, 1).
    C = np.array(rel_errors[:-1]).reshape(-1, 1)
    
    # Use pcolormesh to draw the background
    ax.pcolormesh(
        X, 
        Y, 
        C, 
        cmap=cmap, 
        norm=norm, 
        alpha=alpha, 
        zorder=0, 
        shading='flat'
    )

    ax.set_ylim(lower, upper)


def plot_tolerance_band(
    ax: plt.Axes,
    reference: list[float],
    mode: str,
    tolerance_band: float,
) -> None:
    if tolerance_band <= 0.0:
        return

    if mode == "relative-error":
        ax.axhspan(
            0.0,
            tolerance_band,
            color="#4285f4",
            alpha=0.4,
            linewidth=0.0,
            zorder=1,
        )
        return

    lower = [value * (1.0 - tolerance_band) for value in reference]
    upper = [value * (1.0 + tolerance_band) for value in reference]
    band_low = min(lower + upper)
    band_high = max(lower + upper)
    ax.axhspan(
        band_low,
        band_high,
        color="#4285f4",
        alpha=0.4,
        linewidth=0.0,
        zorder=1,
    )


def center_panel_on_truth(
    ax: plt.Axes,
    truth_values: list[float],
    mode: str,
    y_scale: str,
    tolerance_band: float,
) -> None:
    if tolerance_band <= 0.0 or not truth_values:
        return

    # Use the first truth value as the reference for centering
    truth = truth_values[0]

    if mode == "relative-error":
        center = 0.0
        # Tolerance band is [0, tolerance_band]
        # To make it 50% of the height and center on 0, limits are [-tol, tol]
        low, high = -tolerance_band, tolerance_band
    elif mode == "normalized":
        center = 1.0
        # Tolerance band is [1-tol, 1+tol]. Width is 2*tol.
        # Total height = 4*tol. Limits = 1 +/- 2*tol.
        low, high = 1.0 - 2.0 * tolerance_band, 1.0 + 2.0 * tolerance_band
    else:  # mode == "value"
        center = truth
        # Tolerance band is [truth*(1-tol), truth*(1+tol)]. Width is 2*truth*tol.
        # Total height = 4*truth*tol. Limits = truth +/- 2*truth*tol.
        low, high = truth * (1.0 - 2.0 * tolerance_band), truth * (1.0 + 2.0 * tolerance_band)

    if y_scale == "log" and center > 0:
        # Center in log space: log(high) - log(center) = log(center) - log(low)
        # Range in log space = log(high) - log(low) = 2 * log(high/center)
        # Tolerance band log-width = log(1+tol) - log(1-tol)
        # We want 2 * log(high/center) >= 2 * (log(1+tol) - log(1-tol))
        # log(high/center) >= log((1+tol)/(1-tol))
        ratio = (1.0 + tolerance_band) / (1.0 - tolerance_band)
        low, high = center / ratio, center * ratio

    ax.set_ylim(low, high)


def add_panel(
    ax: plt.Axes,
    x_values: Sequence[float] | Sequence[int],
    estimate: list[float],
    truth: list[float],
    symbol: str,
    panel_label: str,
    mode: str,
    y_scale: str,
    tolerance_band: float,
    tolerance_background: bool,
    tolerance_cmap: Colormap,
    tolerance_norm: Normalize,
    tolerance_background_normalization: str,
    tolerance_background_scale: str,
    tolerance_background_alpha: float,
    highlight_convergence: bool = False,
    center_panel: bool = False,
) -> tuple[plt.Line2D, plt.Line2D]:
    plot_estimate, plot_truth = transform_series(estimate, truth, mode)

    estimate_line, = ax.plot(
        x_values,
        plot_estimate,
        color="#155884",
        linewidth=2.0,
        label="Estimate",
        zorder=3,
    )
    truth_line, = ax.plot(
        x_values,
        plot_truth,
        color="black",
        linewidth=1.8,
        linestyle="--",
        label=reference_label(mode),
        zorder=4,
    )
    configure_y_scale(ax, plot_estimate, plot_truth, rf"${symbol}$", y_scale)
    if tolerance_background:
        plot_tolerance_background(
            ax,
            x_values,
            plot_estimate,
            plot_truth,
            mode,
            tolerance_band,
            tolerance_cmap,
            tolerance_norm,
            tolerance_background_normalization,
            tolerance_background_scale,
            tolerance_background_alpha,
        )
    plot_tolerance_band(ax, plot_truth, mode, tolerance_band)

    if highlight_convergence and tolerance_band > 0.0:
        for x, est, tru in zip(x_values, plot_estimate, plot_truth):
            if tolerance_distance_multiple(est, tru, mode, tolerance_band) <= 1.0:
                ax.axvline(
                    x,
                    color="black",
                    linestyle=":",
                    linewidth=1.5,
                    alpha=0.6,
                    zorder=2,
                )
                break

    if center_panel:
        center_panel_on_truth(ax, plot_truth, mode, y_scale, tolerance_band)

    ax.set_title(rf"${symbol}$")
    ax.text(
        0.99,
        0.98,
        panel_label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_ylabel(y_axis_label(mode, symbol))
    ax.grid(True, alpha=0.3, zorder=2)
    return estimate_line, truth_line


def title_for_mode(mode: str) -> str:
    if mode == "normalized":
        return (
            r"Normalized Estimator Trajectories "
            r"$(\widehat{\sigma}/\sigma, \widehat{\beta}/\beta, "
            r"\widehat{\theta}/\theta)$"
        )
    if mode == "relative-error":
        return (
            r"Relative Error of Estimator Trajectories "
            r"$(\widehat{\sigma}, \widehat{\beta}, \widehat{\theta})$"
        )
    return (
        r"Estimator Trajectories "
        r"$(\widehat{\sigma}, \widehat{\beta}, \widehat{\theta})$"
    )


def main() -> None:
    args = parse_args()
    if args.tolerance_background_max_multiple <= 0.0:
        msg = "--tolerance-background-max-multiple must be positive"
        raise ValueError(msg)
    if not 0.0 <= args.tolerance_background_alpha <= 1.0:
        msg = "--tolerance-background-alpha must be between 0 and 1"
        raise ValueError(msg)

    input_path = resolve_input_path(args.input)
    output_path = resolve_output_path(args.output, input_path)

    raw_trajectory = load_trajectory(input_path)
    min_observation = (
        args.min_observation
        if args.min_observation is not None
        else raw_trajectory["n_obs"][0]
    )
    trajectory = filter_observation_window(
        raw_trajectory,
        min_observation,
        args.max_observation,
    )
    x_values, x_label = choose_x_axis(trajectory, args.x_axis)

    fig, axes = plt.subplots(
        3, 1,
        figsize=(10, 11),
        sharex=True,
    )
    fig.suptitle(
        title_for_mode(args.mode),
        fontsize=15,
        y=0.97,
    )
    tolerance_cmap = plt.get_cmap(args.tolerance_background_cmap)
    
    # Start the colormap at half the tolerance band in percentage
    vmin = args.tolerance_band * 50.0
    
    tolerance_norm = build_tolerance_norm(
        args.tolerance_background_max_multiple,
        args.tolerance_background_scale,
        vmin=vmin,
    )

    estimate_handle = None
    truth_handle = None
    for i, (ax, panel) in enumerate(zip(axes, PARAMETER_PANELS)):
        estimate_handle, truth_handle = add_panel(
            ax,
            x_values,
            trajectory[panel["estimate_key"]],
            trajectory[panel["truth_key"]],
            panel["symbol"],
            panel["panel_label"],
            args.mode,
            args.y_scale,
            args.tolerance_band,
            args.tolerance_background,
            tolerance_cmap,
            tolerance_norm,
            args.tolerance_background_normalization,
            args.tolerance_background_scale,
            args.tolerance_background_alpha,
            highlight_convergence=(i == 0),
            center_panel=(i == 0),
        )

    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[2].set_xlabel(x_label)
    if args.x_axis == "n_obs":
        axes[2].set_xlim(min_observation, args.max_observation)

    legend_handles = [estimate_handle, truth_handle]
    legend_labels = ["Estimate", reference_label(args.mode)]
    if args.tolerance_band > 0.0:
        band_handle = Patch(facecolor="#4285f4", alpha=0.4, edgecolor="#4285f4")
        legend_handles.append(band_handle)
        legend_labels.append(format_tolerance_label(args.tolerance_band))

    show_tolerance_background = args.tolerance_background and args.tolerance_band > 0.0
    if show_tolerance_background:
        colorbar_norm = tolerance_norm
        fig.subplots_adjust(
            left=0.17,
            right=0.96,
            bottom=0.21,
            top=0.92,
            hspace=0.28,
        )
        mappable = ScalarMappable(norm=colorbar_norm, cmap=tolerance_cmap)
        mappable.set_array([])
        colorbar_axis = fig.add_axes([0.30, 0.055, 0.40, 0.018])
        colorbar = fig.colorbar(mappable, cax=colorbar_axis, orientation="horizontal")
        
        colorbar.set_label("Relative error (%)")
        colorbar_ticks = tolerance_colorbar_ticks(
            args.tolerance_background_max_multiple,
            args.tolerance_background_scale,
            vmin=vmin,
        )
        colorbar.set_ticks(colorbar_ticks)
        colorbar.set_ticklabels([f"{t:g}%" for t in colorbar_ticks])
            
        colorbar.ax.xaxis.set_label_position("top")
        colorbar.ax.xaxis.set_ticks_position("bottom")
    else:
        fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))

    fig.legend(
        legend_handles,
        legend_labels,
        loc="lower center",
        ncols=len(legend_handles),
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        edgecolor="0.25",
        bbox_to_anchor=(0.5, 0.105 if show_tolerance_background else 0.02),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
