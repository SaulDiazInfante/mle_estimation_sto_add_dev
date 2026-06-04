#!/usr/bin/env python3
"""ParaView GUI visualization helper for solution snapshot collections.

Load a generated `solution_snapshots.pvd` collection into the ParaView GUI and
create a render view with the configured surface coloring.

Usage:
- Open ParaView GUI
- Tools -> Python Shell
- Run Script... -> select this file

Then call `load_latest_snapshot_collection()` from the Python shell, or edit
`DEFAULT_PRESET` below and run the file directly via the GUI shell.
"""

from pathlib import Path

from paraview.simple import (
    CreateRenderView,
    GetActiveCamera,
    Hide,
    OpenDataFile,
    Render,
    Show,
    WarpByScalar,
)

DEFAULT_INPUT_DIRECTORY = Path("visualization/paraview")
DEFAULT_PRESET = "deterministic-surface"
DEFAULT_IMAGE_SIZE = (1920, 1080)
DEFAULT_BACKGROUND = [1.0, 1.0, 1.0]
DEFAULT_CAMERA_POSITION = [15.0, -20.0, 18.0]
DEFAULT_CAMERA_FOCAL_POINT = [2.5, 2.5, 0.0]
DEFAULT_CAMERA_VIEW_UP = [0.0, 0.0, 1.0]
DEFAULT_WARP_SCALE = 1.0

PRESET_SETTINGS = {
    "stochastic-surface": {
        "scalar": "stochastic",
        "warp_scale": DEFAULT_WARP_SCALE,
        "show_legend": True,
    },
    "deterministic-surface": {
        "scalar": "deterministic",
        "warp_scale": DEFAULT_WARP_SCALE,
        "show_legend": True,
    },
    "flat-deterministic": {
        "scalar": "deterministic",
        "warp_scale": 0.0,
        "show_legend": True,
    },
    "flat-stochastic": {
        "scalar": "stochastic",
        "warp_scale": 0.0,
        "show_legend": True,
    },
}


def find_latest_snapshot_collection(input_dir: Path = DEFAULT_INPUT_DIRECTORY) -> Path:
    """Return the latest solution_snapshots.pvd file under the given directory."""
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    candidates = sorted(input_dir.rglob("solution_snapshots.pvd"))
    if not candidates:
        raise FileNotFoundError(
            f"No solution_snapshots.pvd found under {input_dir}. "
            "Run the export_snapshot_paraview script first."
        )

    return max(candidates, key=lambda path: path.stat().st_mtime)


def configure_camera(camera, position, focal_point, view_up):
    camera.Position = position
    camera.FocalPoint = focal_point
    camera.ViewUp = view_up
    camera.ComputeViewPlaneNormal()


def render_scalar_surface(collection, view, scalar_name, warp_scale, show_legend):
    display = Show(collection, view)
    display.Representation = "Surface"
    display.SetScalarBarVisibility(view, show_legend)
    display.ColorArrayName = ["POINT_DATA", scalar_name]

    if warp_scale == 0.0:
        return collection

    warp = WarpByScalar(Input=collection)
    warp.Scalars = ["POINT_DATA", scalar_name]
    warp.ScaleFactor = warp_scale
    warp_display = Show(warp, view)
    Hide(collection, view)
    warp_display.ColorArrayName = ["POINT_DATA", scalar_name]
    warp_display.SetScalarBarVisibility(view, show_legend)
    return warp


def create_render_view(width, height, background):
    view = CreateRenderView()
    view.ViewSize = [width, height]
    view.Background = background
    view.OrientationAxesVisibility = 1
    return view


def load_snapshot_collection(
    preset: str = DEFAULT_PRESET,
    image_size: tuple[int, int] = DEFAULT_IMAGE_SIZE,
    background: list[float] = DEFAULT_BACKGROUND,
    camera_position: list[float] = DEFAULT_CAMERA_POSITION,
    camera_focal_point: list[float] = DEFAULT_CAMERA_FOCAL_POINT,
    camera_view_up: list[float] = DEFAULT_CAMERA_VIEW_UP,
):
    """Load the latest collection into ParaView GUI and configure the render view."""
    collection_path = find_latest_snapshot_collection()
    print(f"Loading collection: {collection_path}")
    collection = OpenDataFile(str(collection_path))
    if collection is None:
        raise RuntimeError(f"Failed to open input file: {collection_path}")

    settings = PRESET_SETTINGS.get(preset)
    if settings is None:
        raise ValueError(
            f"Unknown preset '{preset}'. Available presets: {', '.join(PRESET_SETTINGS)}"
        )

    width, height = image_size
    view = create_render_view(width, height, background)
    render_scalar_surface(
        collection,
        view,
        settings["scalar"],
        settings["warp_scale"],
        settings["show_legend"],
    )

    configure_camera(GetActiveCamera(view), camera_position,
                     camera_focal_point, camera_view_up)
    Render(view)
    print(f"Created render view for preset '{preset}'.")
    return view


def load_latest_snapshot_collection():
    """Convenience entrypoint for the current default preset."""
    return load_snapshot_collection()


if __name__ == "__main__":
    load_latest_snapshot_collection()
