#!/usr/bin/env python3
"""Plot the time evolution of sigma_hat, beta_hat, and theta_hat."""

from __future__ import annotations

import argparse
import csv
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
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
        default="linear",
        help="Vertical scaling. Use auto to restore sign-aware automatic scaling.",
    )
    parser.add_argument(
        "--tolerance-band",
        type=float,
        default=0.05,
        help="Fractional tolerance band around the reference value; set to 0 to disable.",
    )
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
    print(f"Using symlog on the {title} panel because the series is not strictly positive.")


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


def y_axis_label(mode: str) -> str:
    if mode == "normalized":
        return "Estimate / true value"
    if mode == "relative-error":
        return "Relative error"
    return "Parameter value"


def reference_label(mode: str) -> str:
    if mode == "relative-error":
        return "Zero error"
    if mode == "normalized":
        return "Reference"
    return "True value"


def plot_tolerance_band(
    ax: plt.Axes,
    reference: list[float],
    mode: str,
    tolerance_band: float,
) -> None:
    if tolerance_band <= 0.0:
        return

    if mode == "relative-error":
        ax.axhspan(0.0, tolerance_band, color="#7aa6c2", alpha=0.14, linewidth=0.0)
        return

    lower = [value * (1.0 - tolerance_band) for value in reference]
    upper = [value * (1.0 + tolerance_band) for value in reference]
    band_low = min(lower + upper)
    band_high = max(lower + upper)
    ax.axhspan(band_low, band_high, color="#7aa6c2", alpha=0.14, linewidth=0.0)


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
) -> tuple[plt.Line2D, plt.Line2D]:
    plot_estimate, plot_truth = transform_series(estimate, truth, mode)
    plot_tolerance_band(ax, plot_truth, mode, tolerance_band)

    estimate_line, = ax.plot(
        x_values,
        plot_estimate,
        color="#0b6e4f",
        linewidth=2.0,
        label="Estimate",
    )
    truth_line, = ax.plot(
        x_values,
        plot_truth,
        color="#b22222",
        linewidth=1.8,
        linestyle="--",
        label=reference_label(mode),
    )
    configure_y_scale(ax, plot_estimate, plot_truth, rf"${symbol}$", y_scale)
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
    ax.set_ylabel(y_axis_label(mode))
    ax.grid(True, alpha=0.3)
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
    input_path = resolve_input_path(args.input)
    output_path = resolve_output_path(args.output, input_path)

    trajectory = load_trajectory(input_path)
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

    estimate_handle = None
    truth_handle = None
    for ax, panel in zip(axes, PARAMETER_PANELS):
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
        )

    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[2].set_xlabel(x_label)

    legend_handles = [estimate_handle, truth_handle]
    legend_labels = ["Estimate", reference_label(args.mode)]
    if args.tolerance_band > 0.0:
        band_handle = Patch(facecolor="#7aa6c2", alpha=0.14, edgecolor="none")
        legend_handles.append(band_handle)
        legend_labels.append(f"{args.tolerance_band:.0%} tolerance band")

    fig.legend(
        legend_handles,
        legend_labels,
        loc="lower center",
        ncols=len(legend_handles),
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        edgecolor="0.25",
        bbox_to_anchor=(0.5, 0.02),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
