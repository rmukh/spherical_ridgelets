"""Visualize direct spherical-ridgelet directions in 3D Slicer.

This script converts the custom ``-omd_r`` 24-component NRRD layout into
centered fiber-axis stick glyphs.  It is intended to be run with Slicer's
Python, for example:

  Slicer.exe --python-script tools/slicer/visualize_omd_r.py -- --omd-r out.nrrd

The expected per-voxel component layout is six ``(x, y, z, value)`` entries:

  dir1, -dir1, dir2, -dir2, dir3, -dir3

Only entries 0, 2, and 4 are shown by default because the alternating entries
are antipodes of the same fiber axes.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple


AXIS_COMPONENT_STARTS = (0, 8, 16)


@dataclass
class GlyphRecord:
    i: int
    j: int
    k: int
    center_ras: Tuple[float, float, float]
    direction_ras: Tuple[float, float, float]
    score: float
    axis_id: int


@dataclass
class GlyphSettings:
    slice_axis: str
    slice_index: Optional[int]
    roi_ijk: Optional[Tuple[int, int, int, int, int, int]]
    stride: int
    threshold: float
    max_axes: int
    max_glyphs: int
    glyph_length: Optional[float]
    tube_radius: Optional[float]
    tube_sides: int
    color_by: str
    measurement_frame_mode: str


def _script_args(argv: Sequence[str]) -> List[str]:
    """Accept both normal Python args and Slicer '--python-script x.py -- ...'."""
    args = list(argv[1:])
    if args and args[0] == "--":
        return args[1:]
    return args


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Visualize/export direct ridgelet '-omd_r' directions as stick glyphs."
    )
    parser.add_argument("--omd-r", required=True, help="Input 24-component direct ridgelet maxima NRRD.")
    parser.add_argument("--reference", help="Optional anatomical/dMRI volume to load as background.")
    parser.add_argument("--mask", help="Optional mask volume; voxels <= 0 are skipped.")
    parser.add_argument(
        "--slice-axis",
        default="k",
        choices=("i", "j", "k", "all"),
        help="Voxel axis to slice-filter. Default: k. Use 'all' for whole ROI/volume.",
    )
    parser.add_argument(
        "--slice-index",
        type=int,
        help="Slice index along --slice-axis. Default: middle slice.",
    )
    parser.add_argument(
        "--roi-ijk",
        nargs=6,
        type=int,
        metavar=("IMIN", "IMAX", "JMIN", "JMAX", "KMIN", "KMAX"),
        help="Optional inclusive IJK ROI bounds applied before slice/stride filtering.",
    )
    parser.add_argument("--stride", type=int, default=2, help="Voxel stride in i, j, and k. Default: 2.")
    parser.add_argument("--threshold", type=float, default=0.0, help="Minimum direction score. Default: 0.")
    parser.add_argument("--max-axes", type=int, default=3, choices=(1, 2, 3), help="Max axes per voxel.")
    parser.add_argument("--max-glyphs", type=int, default=100000, help="Safety cap for generated glyphs.")
    parser.add_argument(
        "--glyph-length",
        type=float,
        help="Stick length in physical units. Default: 0.8 * minimum voxel spacing.",
    )
    parser.add_argument(
        "--tube-radius",
        type=float,
        help="Tube radius in physical units. Default: 0.025 * glyph length.",
    )
    parser.add_argument("--tube-sides", type=int, default=8, help="Number of sides for tube glyphs.")
    parser.add_argument("--color-by", default="axis", choices=("axis", "score"), help="Model coloring mode.")
    parser.add_argument("--model-name", default="Direct ridgelet axes", help="Base model name in Slicer.")
    parser.add_argument("--export-vtp", help="Optional VTK XML PolyData export path.")
    parser.add_argument("--export-vtk", help="Optional legacy VTK PolyData export path.")
    parser.add_argument("--export-csv", help="Optional CSV export path.")
    parser.add_argument("--no-display", action="store_true", help="Convert/export without adding model nodes.")
    parser.add_argument(
        "--keep-omd-volume-visible",
        action="store_true",
        help="Keep the loaded 24-component NRRD visible in Slicer's data tree/display.",
    )
    parser.add_argument(
        "--measurement-frame",
        default="auto",
        choices=("auto", "ignore", "require"),
        help="Apply NRRD measurement-frame metadata when found. Default: auto.",
    )
    return parser


def _require_slicer_modules():
    try:
        import numpy as np  # noqa: F401
        import slicer  # noqa: F401
        import vtk  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "This script must be run with 3D Slicer's Python. "
            "Use: Slicer.exe --python-script tools/slicer/visualize_omd_r.py -- --omd-r file.nrrd"
        ) from exc


def _load_volume(path: str, name: Optional[str] = None, vector: bool = False):
    import slicer

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    properties = {}
    if name:
        properties["name"] = name

    if vector:
        try:
            node = slicer.util.loadNodeFromFile(
                path,
                "VectorVolumeFile",
                properties=properties,
                returnNode=True,
            )
            if isinstance(node, tuple):
                success, volume_node = node
                if success:
                    return volume_node
            elif node is not None:
                return node
        except Exception as exc:
            print(f"Warning: VectorVolumeFile load failed, falling back to generic volume loader: {exc}")

    node = slicer.util.loadVolume(path, properties=properties, returnNode=True)
    if isinstance(node, tuple):
        success, volume_node = node
        if not success:
            raise RuntimeError(f"Failed to load volume: {path}")
        return volume_node
    if node is None:
        raise RuntimeError(f"Failed to load volume: {path}")
    return node


def _volume_array_components_last(volume_node):
    import numpy as np
    import slicer

    array = slicer.util.arrayFromVolume(volume_node)
    image_data = volume_node.GetImageData()
    dims_ijk = image_data.GetDimensions()
    expected_kji = (dims_ijk[2], dims_ijk[1], dims_ijk[0])

    if array.ndim != 4:
        raise ValueError(
            f"Expected a 4D array for a 24-component VectorImage, got shape {array.shape}."
        )
    if array.shape[:3] == expected_kji and array.shape[-1] >= 24:
        return array
    if array.shape[1:] == expected_kji and array.shape[0] >= 24:
        return np.moveaxis(array, 0, -1)
    if array.shape[-1] >= 24:
        print(
            f"Warning: unexpected array shape {array.shape}; assuming the last axis stores components."
        )
        return array
    raise ValueError(f"Expected at least 24 components in the -omd_r volume, got shape {array.shape}.")


def _mask_array_kji(mask_node):
    import slicer

    array = slicer.util.arrayFromVolume(mask_node)
    if array.ndim == 4:
        array = array[..., 0]
    if array.ndim != 3:
        raise ValueError(f"Expected a 3D mask volume, got shape {array.shape}.")
    return array


def _vtk_matrix_to_list(matrix) -> List[List[float]]:
    return [[matrix.GetElement(r, c) for c in range(4)] for r in range(4)]


def _ijk_to_ras_matrix(volume_node) -> List[List[float]]:
    import vtk

    matrix = vtk.vtkMatrix4x4()
    volume_node.GetIJKToRASMatrix(matrix)
    return _vtk_matrix_to_list(matrix)


def _point_ijk_to_ras(matrix: Sequence[Sequence[float]], i: int, j: int, k: int):
    x = matrix[0][0] * i + matrix[0][1] * j + matrix[0][2] * k + matrix[0][3]
    y = matrix[1][0] * i + matrix[1][1] * j + matrix[1][2] * k + matrix[1][3]
    z = matrix[2][0] * i + matrix[2][1] * j + matrix[2][2] * k + matrix[2][3]
    return (x, y, z)


def _normalize(vec: Sequence[float]) -> Optional[Tuple[float, float, float]]:
    norm = math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])
    if norm <= 0.0 or not math.isfinite(norm):
        return None
    return (vec[0] / norm, vec[1] / norm, vec[2] / norm)


def _transform_direction(
    direction: Sequence[float], measurement_frame: Optional[Sequence[Sequence[float]]]
) -> Optional[Tuple[float, float, float]]:
    if measurement_frame is None:
        return _normalize(direction)

    transformed = (
        measurement_frame[0][0] * direction[0]
        + measurement_frame[0][1] * direction[1]
        + measurement_frame[0][2] * direction[2],
        measurement_frame[1][0] * direction[0]
        + measurement_frame[1][1] * direction[1]
        + measurement_frame[1][2] * direction[2],
        measurement_frame[2][0] * direction[0]
        + measurement_frame[2][1] * direction[1]
        + measurement_frame[2][2] * direction[2],
    )
    return _normalize(transformed)


def _parse_matrix_from_text(text: str) -> Optional[List[List[float]]]:
    numbers = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    if len(numbers) < 9:
        return None
    values = [float(v) for v in numbers[:9]]
    return [values[0:3], values[3:6], values[6:9]]


def _node_attribute_values(node) -> Iterable[Tuple[str, str]]:
    import vtk

    names = vtk.vtkStringArray()
    if hasattr(node, "GetAttributeNames"):
        node.GetAttributeNames(names)
        for idx in range(names.GetNumberOfValues()):
            key = names.GetValue(idx)
            value = node.GetAttribute(key)
            if value:
                yield key, value


def _measurement_frame_from_node(volume_node) -> Optional[List[List[float]]]:
    known_keys = (
        "NRRD_measurement frame",
        "NRRD_measurement_frame",
        "measurement frame",
        "MeasurementFrame",
        "DWMRI_measurement_frame",
        "DWMRI_measurement-frame",
    )

    nodes_to_check = [volume_node]
    storage_node = volume_node.GetStorageNode() if hasattr(volume_node, "GetStorageNode") else None
    if storage_node is not None:
        nodes_to_check.append(storage_node)

    for node in nodes_to_check:
        for key in known_keys:
            value = node.GetAttribute(key) if hasattr(node, "GetAttribute") else None
            if value:
                matrix = _parse_matrix_from_text(value)
                if matrix is not None:
                    print(f"Using measurement frame from attribute '{key}'.")
                    return matrix
        for key, value in _node_attribute_values(node):
            if "measurement" in key.lower() and "frame" in key.lower():
                matrix = _parse_matrix_from_text(value)
                if matrix is not None:
                    print(f"Using measurement frame from attribute '{key}'.")
                    return matrix
    return None


def _axis_ranges(
    dims_ijk: Tuple[int, int, int], settings: GlyphSettings
) -> Tuple[range, range, range]:
    i_min, i_max = 0, dims_ijk[0] - 1
    j_min, j_max = 0, dims_ijk[1] - 1
    k_min, k_max = 0, dims_ijk[2] - 1

    if settings.roi_ijk is not None:
        ri_min, ri_max, rj_min, rj_max, rk_min, rk_max = settings.roi_ijk
        i_min, i_max = max(i_min, ri_min), min(i_max, ri_max)
        j_min, j_max = max(j_min, rj_min), min(j_max, rj_max)
        k_min, k_max = max(k_min, rk_min), min(k_max, rk_max)

    if i_min > i_max or j_min > j_max or k_min > k_max:
        raise ValueError("The requested ROI does not overlap the -omd_r volume.")

    stride = max(settings.stride, 1)
    if settings.slice_axis != "all":
        if settings.slice_axis == "i":
            idx = settings.slice_index if settings.slice_index is not None else (i_min + i_max) // 2
            i_min = i_max = min(max(idx, i_min), i_max)
        elif settings.slice_axis == "j":
            idx = settings.slice_index if settings.slice_index is not None else (j_min + j_max) // 2
            j_min = j_max = min(max(idx, j_min), j_max)
        elif settings.slice_axis == "k":
            idx = settings.slice_index if settings.slice_index is not None else (k_min + k_max) // 2
            k_min = k_max = min(max(idx, k_min), k_max)

    return (
        range(i_min, i_max + 1, stride),
        range(j_min, j_max + 1, stride),
        range(k_min, k_max + 1, stride),
    )


def collect_glyph_records(volume_node, mask_node, settings: GlyphSettings) -> List[GlyphRecord]:
    array = _volume_array_components_last(volume_node)
    dims_ijk = volume_node.GetImageData().GetDimensions()
    if array.shape[:3] != (dims_ijk[2], dims_ijk[1], dims_ijk[0]):
        raise ValueError(
            f"Array shape {array.shape[:3]} does not match image dimensions {dims_ijk}."
        )

    mask_array = _mask_array_kji(mask_node) if mask_node is not None else None
    if mask_array is not None and mask_array.shape != array.shape[:3]:
        raise ValueError(
            f"Mask shape {mask_array.shape} does not match -omd_r shape {array.shape[:3]}."
        )

    measurement_frame = None
    if settings.measurement_frame_mode != "ignore":
        measurement_frame = _measurement_frame_from_node(volume_node)
        if measurement_frame is None:
            message = (
                "No NRRD measurement-frame metadata found. Assuming stored directions are already "
                "in RAS/physical coordinates."
            )
            if settings.measurement_frame_mode == "require":
                raise ValueError(message)
            print("Warning: " + message)

    ijk_to_ras = _ijk_to_ras_matrix(volume_node)
    i_range, j_range, k_range = _axis_ranges(dims_ijk, settings)

    records: List[GlyphRecord] = []
    for k in k_range:
        for j in j_range:
            for i in i_range:
                if mask_array is not None and mask_array[k, j, i] <= 0:
                    continue

                center = _point_ijk_to_ras(ijk_to_ras, i, j, k)
                for axis_id, start in enumerate(AXIS_COMPONENT_STARTS[: settings.max_axes]):
                    dx, dy, dz, score = array[k, j, i, start : start + 4]
                    score = float(score)
                    if score <= settings.threshold or not math.isfinite(score):
                        continue

                    direction = _transform_direction((float(dx), float(dy), float(dz)), measurement_frame)
                    if direction is None:
                        continue

                    records.append(
                        GlyphRecord(
                            i=i,
                            j=j,
                            k=k,
                            center_ras=center,
                            direction_ras=direction,
                            score=score,
                            axis_id=axis_id,
                        )
                    )
                    if len(records) >= settings.max_glyphs:
                        print(
                            f"Warning: reached --max-glyphs={settings.max_glyphs}; "
                            "increase the limit or use stricter filtering for more glyphs."
                        )
                        return records
    return records


def _record_endpoints(record: GlyphRecord, glyph_length: float):
    half = 0.5 * glyph_length
    cx, cy, cz = record.center_ras
    dx, dy, dz = record.direction_ras
    p0 = (cx - dx * half, cy - dy * half, cz - dz * half)
    p1 = (cx + dx * half, cy + dy * half, cz + dz * half)
    return p0, p1


def records_to_polydata(records: Sequence[GlyphRecord], glyph_length: float, tube_radius: float, tube_sides: int):
    import vtk

    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    score_array = vtk.vtkFloatArray()
    score_array.SetName("score")

    axis_array = vtk.vtkIntArray()
    axis_array.SetName("axis_id")

    direction_array = vtk.vtkFloatArray()
    direction_array.SetName("direction")
    direction_array.SetNumberOfComponents(3)

    ijk_array = vtk.vtkIntArray()
    ijk_array.SetName("ijk")
    ijk_array.SetNumberOfComponents(3)

    for record in records:
        p0, p1 = _record_endpoints(record, glyph_length)
        pid0 = points.InsertNextPoint(p0)
        pid1 = points.InsertNextPoint(p1)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, pid0)
        line.GetPointIds().SetId(1, pid1)
        lines.InsertNextCell(line)

        score_array.InsertNextValue(record.score)
        axis_array.InsertNextValue(record.axis_id)
        direction_array.InsertNextTuple3(*record.direction_ras)
        ijk_array.InsertNextTuple3(record.i, record.j, record.k)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetCellData().AddArray(score_array)
    polydata.GetCellData().AddArray(axis_array)
    polydata.GetCellData().AddArray(direction_array)
    polydata.GetCellData().AddArray(ijk_array)
    polydata.GetCellData().SetActiveScalars("score")

    if tube_radius > 0.0 and records:
        tube = vtk.vtkTubeFilter()
        tube.SetInputData(polydata)
        tube.SetRadius(tube_radius)
        tube.SetNumberOfSides(max(tube_sides, 3))
        tube.CappingOn()
        tube.Update()
        tubed = vtk.vtkPolyData()
        tubed.DeepCopy(tube.GetOutput())
        return tubed

    return polydata


def save_polydata(polydata, export_vtp: Optional[str], export_vtk: Optional[str]) -> None:
    import vtk

    if export_vtp:
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(export_vtp)
        writer.SetInputData(polydata)
        if writer.Write() != 1:
            raise RuntimeError(f"Failed to write VTP: {export_vtp}")
        print(f"Wrote VTP: {export_vtp}")

    if export_vtk:
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(export_vtk)
        writer.SetInputData(polydata)
        if writer.Write() != 1:
            raise RuntimeError(f"Failed to write VTK: {export_vtk}")
        print(f"Wrote VTK: {export_vtk}")


def save_csv(records: Sequence[GlyphRecord], export_csv: str) -> None:
    with open(export_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("i", "j", "k", "x", "y", "z", "axis_id", "dx", "dy", "dz", "score"))
        for record in records:
            writer.writerow(
                (
                    record.i,
                    record.j,
                    record.k,
                    record.center_ras[0],
                    record.center_ras[1],
                    record.center_ras[2],
                    record.axis_id,
                    record.direction_ras[0],
                    record.direction_ras[1],
                    record.direction_ras[2],
                    record.score,
                )
            )
    print(f"Wrote CSV: {export_csv}")


def _add_model_node(polydata, name: str, color: Tuple[float, float, float], opacity: float = 1.0):
    import slicer

    model_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", name)
    model_node.SetAndObservePolyData(polydata)

    display_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode", name + " display")
    display_node.SetColor(color)
    display_node.SetOpacity(opacity)
    display_node.SetVisibility(True)
    model_node.SetAndObserveDisplayNodeID(display_node.GetID())
    return model_node


def display_records(
    records: Sequence[GlyphRecord],
    glyph_length: float,
    tube_radius: float,
    tube_sides: int,
    model_name: str,
    color_by: str,
) -> None:
    if not records:
        print("Warning: no glyphs passed filtering; no model nodes were created.")
        return

    if color_by == "axis":
        axis_colors = (
            (0.95, 0.15, 0.15),
            (0.10, 0.70, 0.25),
            (0.20, 0.35, 1.00),
        )
        for axis_id, color in enumerate(axis_colors):
            axis_records = [record for record in records if record.axis_id == axis_id]
            if not axis_records:
                continue
            polydata = records_to_polydata(axis_records, glyph_length, tube_radius, tube_sides)
            _add_model_node(polydata, f"{model_name} axis {axis_id + 1}", color)
    else:
        polydata = records_to_polydata(records, glyph_length, tube_radius, tube_sides)
        model_node = _add_model_node(polydata, model_name, (1.0, 0.85, 0.15))
        display_node = model_node.GetDisplayNode()
        if display_node is not None:
            scores = [record.score for record in records]
            display_node.SetScalarVisibility(True)
            display_node.SetActiveScalarName("score")
            display_node.SetScalarRange(min(scores), max(scores))


def _hide_volume(volume_node) -> None:
    display_node = volume_node.GetDisplayNode()
    if display_node is not None:
        display_node.SetVisibility(False)


def _show_reference(reference_node) -> None:
    import slicer

    try:
        slicer.util.setSliceViewerLayers(background=reference_node)
    except Exception as exc:
        print(f"Warning: could not set reference volume as slice background: {exc}")


def run(args: argparse.Namespace) -> int:
    _require_slicer_modules()

    settings = GlyphSettings(
        slice_axis=args.slice_axis,
        slice_index=args.slice_index,
        roi_ijk=tuple(args.roi_ijk) if args.roi_ijk else None,
        stride=max(args.stride, 1),
        threshold=args.threshold,
        max_axes=args.max_axes,
        max_glyphs=args.max_glyphs,
        glyph_length=args.glyph_length,
        tube_radius=args.tube_radius,
        tube_sides=args.tube_sides,
        color_by=args.color_by,
        measurement_frame_mode=args.measurement_frame,
    )

    omd_node = _load_volume(args.omd_r, name="direct ridgelet directions", vector=True)
    if not args.keep_omd_volume_visible:
        _hide_volume(omd_node)

    reference_node = _load_volume(args.reference, name="reference") if args.reference else None
    if reference_node is not None:
        _show_reference(reference_node)

    mask_node = _load_volume(args.mask, name="glyph mask") if args.mask else None

    spacing = omd_node.GetSpacing()
    auto_length = 0.8 * min(spacing)
    glyph_length = settings.glyph_length if settings.glyph_length is not None else auto_length
    tube_radius = settings.tube_radius if settings.tube_radius is not None else 0.025 * glyph_length

    records = collect_glyph_records(omd_node, mask_node, settings)
    print(f"Generated {len(records)} glyphs.")

    combined_polydata = records_to_polydata(records, glyph_length, tube_radius, settings.tube_sides)
    save_polydata(combined_polydata, args.export_vtp, args.export_vtk)

    if args.export_csv:
        save_csv(records, args.export_csv)

    if not args.no_display:
        display_records(
            records,
            glyph_length,
            tube_radius,
            settings.tube_sides,
            args.model_name,
            settings.color_by,
        )

    return 0


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(_script_args(argv))
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
