#!/usr/bin/env python3
"""Animate deterministic and stochastic solution snapshots over time."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt

from plot_solution_snapshot_comparison import DEFAULT_PLOT_DIR
from plot_solution_snapshot_comparison import add_velocity_quiver
from plot_solution_snapshot_comparison import centers_to_edges
from plot_solution_snapshot_comparison import choose_color_scale
from plot_solution_snapshot_comparison import compute_velocity_field
from plot_solution_snapshot_comparison import format_time_label
from plot_solution_snapshot_comparison import load_snapshot_data
from plot_solution_snapshot_comparison import read_float_setting
from plot_solution_snapshot_comparison import read_int_setting
from plot_solution_snapshot_comparison import resolve_input_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a time-evolution video comparing deterministic and stochastic snapshots."
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
        help="Output animation path; defaults to MP4 when ffmpeg is available, otherwise GIF.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=8,
        help="Animation frames per second.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="Rendered animation resolution in dots per inch.",
    )
    parser.add_argument(
        "--writer",
        choices=["auto", "ffmpeg", "pillow"],
        default="auto",
        help="Animation writer backend.",
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


def resolve_output_path(argument: str | None, input_path: Path, writer_name: str) -> Path:
    if argument is not None:
        return Path(argument)

    extension = ".mp4" if writer_name == "ffmpeg" else ".gif"
    return DEFAULT_PLOT_DIR / f"{input_path.stem}{extension}"


def choose_writer_name(requested_writer: str, output_argument: str | None) -> str:
    available_writers = set(animation.writers.list())

    if requested_writer == "auto":
        if output_argument is not None and Path(output_argument).suffix.lower() == ".gif":
            if "pillow" in available_writers:
                return "pillow"
            msg = "GIF output requires the pillow animation writer."
            raise RuntimeError(msg)
        if "ffmpeg" in available_writers:
            return "ffmpeg"
        if "pillow" in available_writers:
            return "pillow"
        msg = "No supported animation writer is available. Install ffmpeg or pillow."
        raise RuntimeError(msg)

    if requested_writer not in available_writers:
        msg = f"Requested animation writer '{requested_writer}' is not available."
        raise RuntimeError(msg)

    return requested_writer


def build_writer(writer_name: str, fps: int) -> animation.AbstractMovieWriter:
    if writer_name == "ffmpeg":
        return animation.FFMpegWriter(fps=fps, codec="libx264")
    if writer_name == "pillow":
        return animation.PillowWriter(fps=fps)

    msg = f"Unsupported writer '{writer_name}'."
    raise RuntimeError(msg)


def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)
    writer_name = choose_writer_name(args.writer, args.output)
    output_path = resolve_output_path(args.output, input_path, writer_name)

    if writer_name == "pillow" and output_path.suffix.lower() != ".gif":
        msg = "The pillow writer requires a .gif output path."
        raise ValueError(msg)

    times, x_coordinates, y_coordinates, fields, metadata = load_snapshot_data(input_path)
    x_edges = centers_to_edges(x_coordinates)
    y_edges = centers_to_edges(y_coordinates)
    vmin, vmax, cmap = choose_color_scale(fields)

    velocity_mode_x = read_int_setting(
        args.velocity_mode_x,
        metadata,
        "velocity_mode_x",
        "SARGAZO_VELOCITY_MODE_X",
        1,
    )
    velocity_mode_y = read_int_setting(
        args.velocity_mode_y,
        metadata,
        "velocity_mode_y",
        "SARGAZO_VELOCITY_MODE_Y",
        1,
    )
    length_x = read_float_setting(
        args.length_x,
        metadata,
        "length_x",
        "SARGAZO_LENGTH_X",
        5.0,
    )
    length_y = read_float_setting(
        args.length_y,
        metadata,
        "length_y",
        "SARGAZO_LENGTH_Y",
        5.0,
    )
    quiver_stride = args.quiver_stride
    if quiver_stride is None:
        quiver_stride = max(1, int(max(len(x_coordinates), len(y_coordinates)) / 18))

    x_grid, y_grid, velocity_x, velocity_y = compute_velocity_field(
        x_coordinates,
        y_coordinates,
        velocity_mode_x,
        velocity_mode_y,
        length_x,
        length_y,
    )

    initial_time = times[0]
    deterministic_field = fields[("deterministic", initial_time)]
    stochastic_field = fields[("stochastic", initial_time)]
    extent = (x_edges[0], x_edges[-1], y_edges[0], y_edges[-1])

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.3), sharex=True, sharey=True)
    deterministic_image = axes[0].imshow(
        deterministic_field,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )
    stochastic_image = axes[1].imshow(
        stochastic_field,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )

    for axis, title in zip(axes, ["Deterministic", "Stochastic"]):
        add_velocity_quiver(
            axis,
            x_grid,
            y_grid,
            velocity_x,
            velocity_y,
            quiver_stride,
        )
        axis.set_title(title)
        axis.set_xlabel(r"$x$")
        axis.tick_params(direction="out", length=3.0)
    axes[0].set_ylabel(r"$y$")

    fig.suptitle(
        rf"Time evolution comparison at $t = {format_time_label(initial_time)}$",
        fontsize=15,
    )
    colorbar = fig.colorbar(
        deterministic_image,
        ax=axes,
        fraction=0.034,
        pad=0.035,
    )
    colorbar.set_label("Solution value")

    fig.subplots_adjust(
        left=0.08,
        right=0.88,
        bottom=0.12,
        top=0.88,
        wspace=0.12,
    )

    def update(frame_index: int) -> tuple[plt.AxesImage, plt.AxesImage]:
        time_value = times[frame_index]
        deterministic_image.set_data(fields[("deterministic", time_value)])
        stochastic_image.set_data(fields[("stochastic", time_value)])
        fig.suptitle(
            rf"Time evolution comparison at $t = {format_time_label(time_value)}$",
            fontsize=15,
        )
        return deterministic_image, stochastic_image

    movie = animation.FuncAnimation(
        fig,
        update,
        frames=len(times),
        interval=1000.0 / args.fps,
        blit=False,
        repeat=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    movie.save(output_path, writer=build_writer(writer_name, args.fps), dpi=args.dpi)
    plt.close(fig)

    print(f"Wrote animation to {output_path}")


if __name__ == "__main__":
    main()
