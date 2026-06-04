#!/usr/bin/env python3
"""Export initial/final ParaView surface images for deterministic and stochastic fields.

The script reads a generated `solution_snapshots.pvd` collection and writes two
PNG files:

* `<prefix>_deterministic.png`
* `<prefix>_stochastic.png`

Each image contains the initial and final simulation states for that field.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont, ImageOps

from paraview_solution_snapshot_visualization import (
    DEFAULT_INPUT_DIRECTORY,
    fail,
    find_latest_snapshot_collection,
    load_paraview_api,
)

DEFAULT_IMAGE_SIZE = (2400, 1000)
DEFAULT_CAMERA_POSITION = [7.8, -8.5, 4.8]
DEFAULT_CAMERA_FOCAL_POINT = [2.5, 2.5, 0.25]
DEFAULT_CAMERA_VIEW_UP = [0.0, 0.0, 1.0]
DEFAULT_WARP_SCALE = 1.0
DEFAULT_DETERMINISTIC_WARP_SCALE = 1.0
DEFAULT_STOCHASTIC_WARP_SCALE = 0.8
REPRESENTATION_CHOICES = {
    "surface": "Surface",
    "surface-with-edges": "Surface With Edges",
    "wireframe": "Wireframe",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render deterministic and stochastic initial/final solution "
            "comparisons from a ParaView .pvd collection."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help=(
            "Path to solution_snapshots.pvd. If omitted, the latest collection "
            "under visualization/paraview is used."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="visualization/paraview",
        help="Directory where the two PNG images will be written.",
    )
    parser.add_argument(
        "--output-prefix",
        default="solution_initial_final",
        help="Filename prefix for the deterministic and stochastic PNG files.",
    )
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        default=list(DEFAULT_IMAGE_SIZE),
        help="Final combined PNG resolution as WIDTH HEIGHT.",
    )
    parser.add_argument(
        "--warp-scale",
        type=float,
        default=DEFAULT_WARP_SCALE,
        help="Base scale factor for Warp By Scalar to lift the field into 3-D.",
    )
    parser.add_argument(
        "--deterministic-warp-scale",
        type=float,
        default=None,
        help=(
            "Optional separate warp scale for deterministic panels. "
            "If omitted, --warp-scale is used."
        ),
    )
    parser.add_argument(
        "--stochastic-warp-scale",
        type=float,
        default=None,
        help=(
            "Optional separate warp scale for stochastic panels. "
            "If omitted, --warp-scale is used."
        ),
    )
    parser.add_argument(
        "--representation",
        choices=sorted(REPRESENTATION_CHOICES),
        default="surface-with-edges",
        help="ParaView surface representation used for the rendered fields.",
    )
    parser.add_argument(
        "--camera-position",
        type=float,
        nargs=3,
        default=DEFAULT_CAMERA_POSITION,
        help="Camera position in world coordinates.",
    )
    parser.add_argument(
        "--camera-focal-point",
        type=float,
        nargs=3,
        default=DEFAULT_CAMERA_FOCAL_POINT,
        help="Camera focal point in world coordinates.",
    )
    parser.add_argument(
        "--camera-view-up",
        type=float,
        nargs=3,
        default=DEFAULT_CAMERA_VIEW_UP,
        help="Camera view-up vector.",
    )
    parser.add_argument(
        "--background",
        type=float,
        nargs=3,
        default=[1.0, 1.0, 1.0],
        help="Background color as RGB values between 0 and 1.",
    )
    parser.add_argument(
        "--hide-legend",
        action="store_true",
        help="Hide the scalar-bar legend in the rendered panels.",
    )
    return parser.parse_args()


def resolve_input_path(argument: str | None) -> Path:
    if argument:
        input_path = Path(argument)
    else:
        input_path = find_latest_snapshot_collection(DEFAULT_INPUT_DIRECTORY)
        if input_path is not None:
            print(f"Using latest input file: {input_path}")

    if input_path is None or not input_path.exists():
        fail(
            "Input file does not exist. Provide --input or generate a snapshot "
            "collection under visualization/paraview with solution_snapshots.pvd."
        )

    return input_path


def get_time_steps(scene: Any) -> list[float]:
    scene.UpdateAnimationUsingDataTimeSteps()
    timekeeper = scene.TimeKeeper
    if hasattr(timekeeper, "TimeSteps"):
        values = list(timekeeper.TimeSteps)
    elif hasattr(timekeeper, "TimestepValues"):
        values = list(timekeeper.TimestepValues)
    else:
        values = []

    if not values:
        fail("Could not read time steps from the ParaView animation scene.")

    return values


def configure_camera(api: Any, view: Any, args: argparse.Namespace) -> None:
    view.CameraPosition = args.camera_position
    view.CameraFocalPoint = args.camera_focal_point
    view.CameraViewUp = args.camera_view_up
    view.CameraParallelProjection = 0
    view.CameraViewAngle = 32.0

    if hasattr(view, "GetActiveCamera"):
        camera = view.GetActiveCamera()
    else:
        camera = api.GetActiveCamera()
    camera.Position = args.camera_position
    camera.FocalPoint = args.camera_focal_point
    camera.ViewUp = args.camera_view_up
    camera.ComputeViewPlaneNormal()
    view.CenterOfRotation = args.camera_focal_point


def set_display_style(
    display: Any,
    view: Any,
    scalar_name: str,
    representation: str,
    show_legend: bool,
) -> None:
    display.Representation = representation
    display.ColorArrayName = ["POINT_DATA", scalar_name]
    display.SetScalarBarVisibility(view, show_legend)

    if hasattr(display, "RescaleTransferFunctionToDataRange"):
        display.RescaleTransferFunctionToDataRange(True, False)

    for attr_name, value in [
        ("Ambient", 0.18),
        ("Diffuse", 0.82),
        ("Specular", 0.22),
        ("SpecularPower", 18.0),
    ]:
        try:
            setattr(display, attr_name, value)
        except Exception:
            pass


def render_scalar_surface(
    api: Any,
    collection: Any,
    view: Any,
    scalar_name: str,
    warp_scale: float,
    representation: str,
    show_legend: bool,
) -> Any:
    base_display = api.Show(collection, view)
    set_display_style(
        base_display, view, scalar_name, representation, show_legend
    )

    if warp_scale == 0.0:
        return collection

    warp = api.WarpByScalar(Input=collection)
    warp.Scalars = ["POINT_DATA", scalar_name]
    warp.ScaleFactor = warp_scale
    warp_display = api.Show(warp, view)
    api.Hide(collection, view)
    set_display_style(
        warp_display, view, scalar_name, representation, show_legend
    )
    return warp


def save_timepoint_panel(
    api: Any,
    collection: Any,
    scalar_name: str,
    time_value: float,
    output_path: Path,
    width: int,
    height: int,
    args: argparse.Namespace,
) -> None:
    view = api.CreateRenderView()
    view.ViewSize = [width, height]
    view.Background = args.background
    view.Background2 = args.background
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 1
    view.ViewTime = time_value

    render_scalar_surface(
        api,
        collection,
        view,
        scalar_name,
        get_warp_scale(scalar_name, args),
        REPRESENTATION_CHOICES[args.representation],
        not args.hide_legend,
    )
    configure_camera(api, view, args)
    api.Render(view)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    api.SaveScreenshot(
        str(output_path),
        view,
        ImageResolution=[width, height],
        TransparentBackground=0,
    )


def load_font(size: int) -> ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def format_time(time_value: float) -> str:
    return f"{time_value:.6g}"


def paste_contained(
    canvas: Image.Image,
    image: Image.Image,
    box: tuple[int, int, int, int],
) -> None:
    left, top, right, bottom = box
    target_size = (right - left, bottom - top)
    contained = ImageOps.contain(
        image, target_size, method=Image.Resampling.LANCZOS)
    x = left + (target_size[0] - contained.width) // 2
    y = top + (target_size[1] - contained.height) // 2
    canvas.paste(contained, (x, y))


def compose_comparison_image(
    initial_panel: Path,
    final_panel: Path,
    output_path: Path,
    field_label: str,
    initial_time: float,
    final_time: float,
    size: tuple[int, int],
) -> None:
    width, height = size
    margin = max(24, width // 70)
    gutter = max(18, width // 100)
    title_height = max(58, height // 13)
    label_height = max(42, height // 20)
    panel_width = (width - 2 * margin - gutter) // 2
    panel_height = height - 2 * margin - title_height - label_height

    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(max(24, height // 34))
    label_font = load_font(max(18, height // 48))

    title = f"{field_label} solution: initial and final states"
    draw.text((margin, margin), title, fill=(20, 20, 20), font=title_font)

    left_x = margin
    right_x = margin + panel_width + gutter
    label_y = margin + title_height
    panel_y = label_y + label_height
    left_box = (left_x, panel_y, left_x + panel_width, panel_y + panel_height)
    right_box = (right_x, panel_y, right_x +
                 panel_width, panel_y + panel_height)

    draw.text(
        (left_x, label_y),
        f"Initial time  t = {format_time(initial_time)}",
        fill=(35, 35, 35),
        font=label_font,
    )
    draw.text(
        (right_x, label_y),
        f"Final time  t = {format_time(final_time)}",
        fill=(35, 35, 35),
        font=label_font,
    )

    with Image.open(initial_panel).convert("RGB") as initial_image:
        paste_contained(canvas, initial_image, left_box)
    with Image.open(final_panel).convert("RGB") as final_image:
        paste_contained(canvas, final_image, right_box)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_path)
    print(f"Saved {field_label.lower()} initial/final image: {output_path}")


def get_warp_scale(scalar_name: str, args: argparse.Namespace) -> float:
    if scalar_name == "deterministic":
        return (
            args.deterministic_warp_scale
            if args.deterministic_warp_scale is not None
            else args.warp_scale
        )
    if scalar_name == "stochastic":
        return (
            args.stochastic_warp_scale
            if args.stochastic_warp_scale is not None
            else args.warp_scale
        )
    return args.warp_scale


def export_field_comparison(
    api: Any,
    collection: Any,
    scalar_name: str,
    output_path: Path,
    initial_time: float,
    final_time: float,
    args: argparse.Namespace,
) -> None:
    width, height = args.size
    panel_width = max(640, width // 2)
    panel_height = max(480, height)

    with tempfile.TemporaryDirectory(prefix="paraview_initial_final_") as temp_name:
        temp_dir = Path(temp_name)
        initial_panel = temp_dir / f"{scalar_name}_initial.png"
        final_panel = temp_dir / f"{scalar_name}_final.png"

        save_timepoint_panel(
            api,
            collection,
            scalar_name,
            initial_time,
            initial_panel,
            panel_width,
            panel_height,
            args,
        )
        save_timepoint_panel(
            api,
            collection,
            scalar_name,
            final_time,
            final_panel,
            panel_width,
            panel_height,
            args,
        )
        compose_comparison_image(
            initial_panel,
            final_panel,
            output_path,
            scalar_name.capitalize(),
            initial_time,
            final_time,
            tuple(args.size),
        )


def main() -> None:
    args = parse_args()
    input_path = resolve_input_path(args.input)
    output_dir = Path(args.output_dir)

    api = load_paraview_api()
    collection = api.OpenDataFile(str(input_path))
    if collection is None:
        fail(f"Could not open input file: {input_path}")

    scene = api.GetAnimationScene()
    time_steps = get_time_steps(scene)
    initial_time = time_steps[0]
    final_time = time_steps[-1]

    outputs = {
        "deterministic": output_dir / f"{args.output_prefix}_deterministic.png",
        "stochastic": output_dir / f"{args.output_prefix}_stochastic.png",
    }
    for scalar_name, output_path in outputs.items():
        export_field_comparison(
            api,
            collection,
            scalar_name,
            output_path,
            initial_time,
            final_time,
            args,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        sys.exit(130)
