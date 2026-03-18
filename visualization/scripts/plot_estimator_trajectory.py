#!/usr/bin/env python3
"""Plot the time evolution of sigma_hat, beta_hat, and theta_hat."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

DEFAULT_INPUT_DIR = Path("data/output")
DEFAULT_PLOT_DIR = Path("visualization/plots")
ESTIMATOR_SUFFIX = "_estimator_trajectory.csv"


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
    title: str
) -> None:
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


def add_panel(
    ax: plt.Axes,
    time: list[float],
    estimate: list[float],
    truth: list[float],
    title: str
) -> None:
    ax.plot(
        time,
        estimate,
        color="#0b6e4f",
        linewidth=2.0,
        label=f"{title} estimate"
    )
    ax.plot(
        time,
        truth,
        color="#b22222",
        linewidth=1.8,
        linestyle="--",
        label=f"{title} true"
    )
    configure_y_scale(ax, estimate, truth, title)
    ax.set_title(title)
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")


def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)
    output_path = resolve_output_path(args.output, input_path)

    trajectory = load_trajectory(input_path)

    fig, axes = plt.subplots(
        3, 1,
        figsize=(10, 11),
        sharex=True,
        constrained_layout=True
    )
    fig.suptitle("Asymptotic Evolution of the Estimators", fontsize=15)

    add_panel(
        axes[0],
        trajectory["time"],
        trajectory["sigma_hat"],
        trajectory["sigma_true"],
        "sigma"
    )
    add_panel(
        axes[1],
        trajectory["time"],
        trajectory["beta_hat"],
        trajectory["beta_true"],
        "beta"
    )
    add_panel(
        axes[2],
        trajectory["time"],
        trajectory["theta_hat"],
        trajectory["theta_true"],
        "theta"
    )

    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[2].set_xlabel("Time")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
