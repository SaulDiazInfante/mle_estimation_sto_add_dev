#!/usr/bin/env python3
"""Export ParaView screenshots and animations from a .pvd snapshot collection.

This script is intended for use with the generated file
`visualization/paraview/<timestamp>_paraview/solution_snapshots.pvd`.
It supports three presets, static screenshot export, and an optional
frame-number overlay for video output.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple, Optional

DEFAULT_INPUT_DIRECTORY = Path("visualization/paraview")
DEFAULT_IMAGE_SIZE = (1920, 1080)
PRESET_CHOICES = (
    "stochastic-surface",
    "deterministic-surface",
    "side-by-side-comparison",
)


class ParaViewAPI(NamedTuple):
    """Container for imported ParaView Python API entrypoints."""

    CreateRenderView: Any
    GetActiveCamera: Any
    GetAnimationScene: Any
    Hide: Any
    OpenDataFile: Any
    Render: Any
    SaveAnimation: Any
    SaveScreenshot: Any
    Show: Any
    WarpByScalar: Any


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Load a ParaView .pvd time collection and export a static screenshot "
            "and/or a time animation video."
        )
    )
    parser.add_argument(
        "--input",
        default=None,
        help=(
            "Path to the input ParaView .pvd collection. If omitted, the script "
            "will use the latest solution_snapshots.pvd found beneath "
            "visualization/paraview."
        ),
    )
    parser.add_argument(
        "--output-image",
        default=None,
        help="Output screenshot image path (png or jpg).",
    )
    parser.add_argument(
        "--output-initial-image",
        default=None,
        help=(
            "Output PNG path for an initial-condition comparison at time zero. "
            "If used with side-by-side mode, the script will generate two files "
            "with deterministic and stochastic suffixes."
        ),
    )
    parser.add_argument(
        "--output-final-image",
        default=None,
        help=(
            "Output PNG path for a final-frame comparison at the last data timestep. "
            "If used with side-by-side mode, the script will generate two files "
            "with deterministic and stochastic suffixes."
        ),
    )
    parser.add_argument(
        "--output-video",
        default=None,
        help="Output animation path (mp4, avi, or gif).",
    )
    parser.add_argument(
        "--preset",
        choices=PRESET_CHOICES,
        default=None,
        help=(
            "Visualization preset. 'stochastic-surface' and "
            "'deterministic-surface' export a single scalar surface. "
            "'side-by-side-comparison' exports separate deterministic and "
            "stochastic surface outputs."
        ),
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=8,
        help="Frames per second for the output animation.",
    )
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        default=list(DEFAULT_IMAGE_SIZE),
        help="Output image/video resolution as WIDTH HEIGHT.",
    )
    parser.add_argument(
        "--scalar",
        default="stochastic",
        help="Point field used for color mapping; 'stochastic' or 'deterministic'.",
    )
    parser.add_argument(
        "--warp-scale",
        type=float,
        default=1.0,
        help=(
            "Scale factor for the Warp By Scalar filter. Use 0 to disable warp "
            "and keep the original plane geometry."
        ),
    )
    parser.add_argument(
        "--camera-position",
        type=float,
        nargs=3,
        default=[15.0, -20.0, 18.0],
        help="Camera position in world coordinates.",
    )
    parser.add_argument(
        "--camera-focal-point",
        type=float,
        nargs=3,
        default=[2.5, 2.5, 0.0],
        help="Camera focal point in world coordinates.",
    )
    parser.add_argument(
        "--camera-view-up",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 1.0],
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
        "--show-legend",
        action="store_true",
        help="Show the scalar bar legend in the screenshot or animation.",
    )
    parser.add_argument(
        "--frame-label",
        action="store_true",
        help=(
            "Overlay frame numbering in the top right of the exported video "
            "using ffmpeg after the movie has been generated."
        ),
    )
    return parser.parse_args()


def find_latest_snapshot_collection(input_dir: Path) -> Optional[Path]:
    """Return the latest solution_snapshots.pvd file found under the input directory."""
    if not input_dir.exists():
        return None

    candidates = sorted(input_dir.rglob("solution_snapshots.pvd"))
    if not candidates:
        return None

    return max(candidates, key=lambda path: path.stat().st_mtime)


def apply_preset(args: argparse.Namespace) -> None:
    """Apply named presets to the parsed arguments."""
    if args.preset == "stochastic-surface":
        args.scalar = "stochastic"
        args.warp_scale = 1.0
        args.show_legend = True
    elif args.preset == "deterministic-surface":
        args.scalar = "deterministic"
        args.warp_scale = 1.0
        args.show_legend = True
    elif args.preset == "side-by-side-comparison":
        args.show_legend = True


def fail(message: str) -> None:
    """Print an error message and exit with a non-zero status."""
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def load_paraview_api() -> ParaViewAPI:
    """Import ParaView Python entrypoints at runtime.

    This script is designed to be executed with pvpython or pvbatch.
    """
    try:
        from paraview.simple import (
            CreateRenderView,
            GetActiveCamera,
            GetAnimationScene,
            Hide,
            OpenDataFile,
            Render,
            SaveAnimation,
            SaveScreenshot,
            Show,
            WarpByScalar,
        )
    except ImportError:
        fail(
            "ParaView Python modules are not available. "
            "Run this script with pvpython or pvbatch from a ParaView installation."
        )

    return ParaViewAPI(
        CreateRenderView=CreateRenderView,
        GetActiveCamera=GetActiveCamera,
        GetAnimationScene=GetAnimationScene,
        Hide=Hide,
        OpenDataFile=OpenDataFile,
        Render=Render,
        SaveAnimation=SaveAnimation,
        SaveScreenshot=SaveScreenshot,
        Show=Show,
        WarpByScalar=WarpByScalar,
    )


def ffmpeg_available() -> bool:
    """Return True when ffmpeg is installed and available on PATH."""
    return shutil.which("ffmpeg") is not None


def overlay_frame_numbers(
    input_video: Path,
    output_video: Path,
    total_frames: int,
    fontsize: int = 24,
    margin: int = 16,
) -> None:
    """Overlay a frame count label on the generated video using ffmpeg."""
    if not ffmpeg_available():
        fail("ffmpeg is required to overlay frame numbers on the final video.")

    fontfile = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
    font_option = f":fontfile={fontfile}" if Path(fontfile).exists() else ""
    drawtext = (
        "drawtext="
        f"text='Frame %{{n+1}}/{total_frames}':"
        f"x=w-tw-{margin}:y={margin}:fontsize={fontsize}:"
        "fontcolor=white:box=1:boxcolor=black@0.5"
        f"{font_option}"
    )
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-i",
            str(input_video),
            "-vf",
            drawtext,
            "-c:a",
            "copy",
            str(output_video),
        ],
        check=True,
    )


def render_scalar_surface(
    api: ParaViewAPI,
    collection: Any,
    view: Any,
    scalar_name: str,
    warp_scale: float,
    show_legend: bool,
) -> Any:
    """Configure a ParaView view to render a scalar surface."""
    display = api.Show(collection, view)
    display.Representation = "Surface"
    display.SetScalarBarVisibility(view, show_legend)
    display.ColorArrayName = ["POINT_DATA", scalar_name]

    if warp_scale == 0.0:
        return collection

    warp = api.WarpByScalar(Input=collection)
    warp.Scalars = ["POINT_DATA", scalar_name]
    warp.ScaleFactor = warp_scale
    warp_display = api.Show(warp, view)
    api.Hide(collection, view)
    warp_display.ColorArrayName = ["POINT_DATA", scalar_name]
    warp_display.SetScalarBarVisibility(view, show_legend)
    return warp


def save_screenshot(
    api: ParaViewAPI,
    view: Any,
    output_path: Path,
    width: int,
    height: int,
) -> None:
    """Save a rendered view to a screenshot file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    api.SaveScreenshot(
        str(output_path),
        view,
        ImageResolution=[width, height],
        TransparentBackground=0,
    )
    print(f"Saved screenshot: {output_path}")


def _ffmpeg_encode_sequence(
    image_pattern: str,
    output_path: Path,
    fps: int,
    work_dir: Path,
) -> None:
    if not ffmpeg_available():
        fail(
            "ffmpeg is required to generate a video file when ParaView cannot "
            "write the requested movie format directly."
        )

    output_ext = output_path.suffix.lower()
    if output_ext == ".mp4":
        cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            image_pattern,
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            str(output_path),
        ]
    elif output_ext == ".avi":
        cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            image_pattern,
            "-c:v",
            "mpeg4",
            "-qscale:v",
            "2",
            str(output_path),
        ]
    elif output_ext == ".ogv":
        cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            image_pattern,
            "-c:v",
            "libtheora",
            "-qscale:v",
            "5",
            str(output_path),
        ]
    elif output_ext == ".gif":
        palette = work_dir / "palette.png"
        subprocess.run(
            [
                "ffmpeg",
                "-y",
                "-framerate",
                str(fps),
                "-i",
                image_pattern,
                "-vf",
                "palettegen",
                str(palette),
            ],
            check=True,
        )
        subprocess.run(
            [
                "ffmpeg",
                "-y",
                "-framerate",
                str(fps),
                "-i",
                image_pattern,
                "-i",
                str(palette),
                "-lavfi",
                "paletteuse",
                str(output_path),
            ],
            check=True,
        )
        return
    else:
        fail(
            f"Unsupported video fallback format '{output_ext}'. "
            "Supported formats are .mp4, .avi, .ogv, and .gif."
        )

    subprocess.run(cmd, check=True)


def _save_animation_to_image_sequence(
    api: ParaViewAPI,
    view: Any,
    work_dir: Path,
    width: int,
    height: int,
    fps: int,
) -> str:
    image_template = str(work_dir / "frame_{0:05d}.png")
    api.SaveAnimation(
        image_template,
        view,
        ImageResolution=[width, height],
        FrameRate=fps,
        SuffixFormat="{0:05d}",
    )

    generated_frames = sorted(work_dir.glob("frame_*.png"))
    if not generated_frames:
        fail(
            "ParaView failed to generate an image sequence for the fallback animation. "
            "Check whether images are being written to the temporary directory."
        )

    # Convert an example filename like 'frame_00001.png' into a
    # printf-style pattern 'frame_%05d.png' that ffmpeg accepts.
    example_name = generated_frames[0].name
    m = re.match(r"^(.*?)(\d+)(\.[^.]+)?$", example_name)
    if not m:
        pattern_printf = example_name
    else:
        prefix, digits, suffix = m.group(1), m.group(2), m.group(3) or ""
        width_digits = len(digits)
        pattern_printf = f"{prefix}%0{width_digits}d{suffix}"

    return str(work_dir / pattern_printf)


def _prefer_ffmpeg_fallback(output_extension: str) -> bool:
    return output_extension == ".mp4"


def save_animation(
    api: ParaViewAPI,
    view: Any,
    output_path: Path,
    width: int,
    height: int,
    fps: int,
    frame_label: bool,
) -> None:
    """Save a rendered view to an animation file and optionally overlay frame numbers."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="paraview_anim_"))
    start_time = time.monotonic()
    print(f"Starting animation export: {output_path}")
    try:
        if frame_label:
            temp_path = temp_dir / output_path.name
        else:
            temp_path = output_path

        scene = api.GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()

        try:
            if _prefer_ffmpeg_fallback(output_path.suffix.lower()):
                raise RuntimeError(
                    "Direct ParaView movie export is disabled for this format.")

            api.SaveAnimation(
                str(temp_path),
                view,
                ImageResolution=[width, height],
                FrameRate=fps,
            )
            if not temp_path.exists():
                raise RuntimeError(
                    "ParaView did not create the requested output file.")
        except Exception:
            print(
                "INFO: Using image sequence + ffmpeg fallback for movie export.",
                file=sys.stderr,
            )
            image_pattern = _save_animation_to_image_sequence(
                api, view, temp_dir, width, height, fps
            )
            _ffmpeg_encode_sequence(image_pattern, temp_path, fps, temp_dir)

        if frame_label:
            timekeeper = scene.TimeKeeper
            if hasattr(timekeeper, "TimeSteps"):
                total_frames = len(timekeeper.TimeSteps)
            elif hasattr(timekeeper, "TimestepValues"):
                total_frames = len(timekeeper.TimestepValues)
            else:
                total_frames = 0

            if total_frames <= 0:
                fail(
                    "Unable to determine the number of animation frames from the ParaView time keeper. "
                    "Frame-label overlay is not available for this build of ParaView."
                )

            overlay_frame_numbers(temp_path, output_path, total_frames)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    elapsed = time.monotonic() - start_time
    print(f"Saved animation: {output_path} ({elapsed:.1f} seconds)")


def save_timepoint_screenshot(
    api: ParaViewAPI,
    scene: Any,
    view: Any,
    output_path: Path,
    width: int,
    height: int,
    time_point: str,
) -> None:
    """Save a screenshot for a specific animation timepoint."""
    scene.UpdateAnimationUsingDataTimeSteps()

    if time_point == "initial":
        scene.GoToFirst()
    elif time_point == "final":
        scene.GoToLast()
    else:
        raise ValueError("time_point must be 'initial' or 'final'")

    api.Render(view)
    save_screenshot(api, view, output_path, width, height)


def configure_camera(api: ParaViewAPI, args: argparse.Namespace) -> None:
    """Configure the active camera for the current view."""
    camera = api.GetActiveCamera()
    camera.Position = args.camera_position
    camera.FocalPoint = args.camera_focal_point
    camera.ViewUp = args.camera_view_up
    camera.ComputeViewPlaneNormal()


def create_side_by_side_paths(
    output_image: Optional[Path], output_video: Optional[Path]
) -> tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path]]:
    """Create distinct deterministic and stochastic output paths for side-by-side exports."""
    if output_image is None and output_video is None:
        return None, None, None, None

    deterministic_image = (
        output_image.with_name(
            f"{output_image.stem}_deterministic{output_image.suffix}")
        if output_image is not None
        else None
    )
    stochastic_image = (
        output_image.with_name(
            f"{output_image.stem}_stochastic{output_image.suffix}")
        if output_image is not None
        else None
    )
    deterministic_video = (
        output_video.with_name(
            f"{output_video.stem}_deterministic{output_video.suffix}")
        if output_video is not None
        else None
    )
    stochastic_video = (
        output_video.with_name(
            f"{output_video.stem}_stochastic{output_video.suffix}")
        if output_video is not None
        else None
    )
    return deterministic_image, stochastic_image, deterministic_video, stochastic_video


def create_side_by_side_timepoint_paths(
    output_path: Optional[Path], time_label: str
) -> tuple[Optional[Path], Optional[Path]]:
    """Create deterministic and stochastic output paths for a specific timepoint."""
    if output_path is None:
        return None, None

    return (
        output_path.with_name(
            f"{output_path.stem}_deterministic_{time_label}{output_path.suffix}"
        ),
        output_path.with_name(
            f"{output_path.stem}_stochastic_{time_label}{output_path.suffix}"
        ),
    )


def main() -> None:
    args = parse_args()
    apply_preset(args)

    if (
        args.output_image is None
        and args.output_initial_image is None
        and args.output_final_image is None
        and args.output_video is None
    ):
        fail(
            "At least one of --output-image, --output-initial-image, "
            "--output-final-image, or --output-video must be provided."
        )

    if args.input:
        input_path = Path(args.input)
    else:
        input_path = find_latest_snapshot_collection(DEFAULT_INPUT_DIRECTORY)
        if input_path is not None:
            print(f"Using latest input file: {input_path}")

    if input_path is None or not input_path.exists():
        fail(
            "Input file does not exist. Provide --input or generate a snapshot collection "
            "under visualization/paraview with solution_snapshots.pvd."
        )

    output_image = Path(args.output_image) if args.output_image else None
    output_initial_image = (
        Path(args.output_initial_image) if args.output_initial_image else None
    )
    output_final_image = (
        Path(args.output_final_image) if args.output_final_image else None
    )
    output_video = Path(args.output_video) if args.output_video else None
    width, height = args.size
    api = load_paraview_api()

    collection = api.OpenDataFile(str(input_path))
    if collection is None:
        fail(f"Could not open input file: {input_path}")

    if args.preset == "side-by-side-comparison":
        deterministic_image, stochastic_image, deterministic_video, stochastic_video = (
            create_side_by_side_paths(output_image, output_video)
        )
        deterministic_initial_image, stochastic_initial_image = (
            create_side_by_side_timepoint_paths(
                output_initial_image, "initial")
        )
        deterministic_final_image, stochastic_final_image = (
            create_side_by_side_timepoint_paths(output_final_image, "final")
        )
        for scalar_name, image_path, initial_path, final_path, video_path in [
            (
                "deterministic",
                deterministic_image,
                deterministic_initial_image,
                deterministic_final_image,
                deterministic_video,
            ),
            (
                "stochastic",
                stochastic_image,
                stochastic_initial_image,
                stochastic_final_image,
                stochastic_video,
            ),
        ]:
            view = api.CreateRenderView()
            view.ViewSize = [width, height]
            view.Background = args.background
            view.OrientationAxesVisibility = 1
            view.CenterOfRotation = args.camera_focal_point
            render_scalar_surface(
                api,
                collection,
                view,
                scalar_name,
                args.warp_scale,
                args.show_legend,
            )
            configure_camera(api, args)
            api.Render(view)
            if image_path is not None:
                save_screenshot(api, view, image_path, width, height)
            scene = api.GetAnimationScene()
            if initial_path is not None:
                save_timepoint_screenshot(
                    api,
                    scene,
                    view,
                    initial_path,
                    width,
                    height,
                    "initial",
                )
            if final_path is not None:
                save_timepoint_screenshot(
                    api,
                    scene,
                    view,
                    final_path,
                    width,
                    height,
                    "final",
                )
            if video_path is not None:
                save_animation(api, view, video_path, width,
                               height, args.fps, args.frame_label)
        return

    view = api.CreateRenderView()
    view.ViewSize = [width, height]
    view.Background = args.background
    view.OrientationAxesVisibility = 1
    view.CenterOfRotation = args.camera_focal_point
    render_scalar_surface(
        api,
        collection,
        view,
        args.scalar,
        args.warp_scale,
        args.show_legend,
    )
    configure_camera(api, args)
    api.Render(view)

    if output_image is not None:
        save_screenshot(api, view, output_image, width, height)
    scene = api.GetAnimationScene()
    if output_initial_image is not None:
        save_timepoint_screenshot(
            api,
            scene,
            view,
            output_initial_image,
            width,
            height,
            "initial",
        )
    if output_final_image is not None:
        save_timepoint_screenshot(
            api,
            scene,
            view,
            output_final_image,
            width,
            height,
            "final",
        )
    if output_video is not None:
        save_animation(api, view, output_video, width,
                       height, args.fps, args.frame_label)


if __name__ == "__main__":
    main()
