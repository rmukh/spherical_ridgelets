"""Core conversion and visualization functions for RidgeletDirections."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import slicer
import vtk


AXIS_COMPONENT_STARTS = (0, 8, 16)
GENERATED_MODEL_ATTRIBUTE = "RidgeletDirections.GeneratedModel"


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
    slice_axis: str = "k"
    slice_index: Optional[int] = None
    roi_ijk: Optional[Tuple[int, int, int, int, int, int]] = None
    stride: int = 2
    threshold: float = 0.0
    max_axes: int = 3
    max_glyphs: int = 100000
    glyph_length: Optional[float] = None
    tube_radius: Optional[float] = None
    tube_sides: int = 8
    color_by: str = "axis"
    measurement_frame_mode: str = "auto"
    clear_previous_models: bool = True


@dataclass
class VisualizationResult:
    glyph_count: int
    model_nodes: List
    glyph_length: float
    tube_radius: float
    warnings: List[str]


def volume_array_components_last(volume_node):
    import numpy as np

    array = slicer.util.arrayFromVolume(volume_node)
    dims_ijk = volume_node.GetImageData().GetDimensions()
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
        return array
    raise ValueError(f"Expected at least 24 components in the -omd_r volume, got shape {array.shape}.")


def mask_array_kji(mask_node):
    array = slicer.util.arrayFromVolume(mask_node)
    if array.ndim == 4:
        array = array[..., 0]
    if array.ndim != 3:
        raise ValueError(f"Expected a 3D mask volume, got shape {array.shape}.")
    return array


def validate_input_volume(volume_node) -> None:
    if volume_node is None or volume_node.GetImageData() is None:
        raise ValueError("Select a direct ridgelet directions volume.")
    components = volume_node.GetImageData().GetNumberOfScalarComponents()
    if components < 24:
        raise ValueError(
            f"The selected volume has {components} components; an -omd_r volume needs at least 24."
        )


def _vtk_matrix_to_list(matrix) -> List[List[float]]:
    return [[matrix.GetElement(r, c) for c in range(4)] for r in range(4)]


def _ijk_to_ras_matrix(volume_node) -> List[List[float]]:
    matrix = vtk.vtkMatrix4x4()
    volume_node.GetIJKToRASMatrix(matrix)
    return _vtk_matrix_to_list(matrix)


def _point_ijk_to_ras(matrix: Sequence[Sequence[float]], i: int, j: int, k: int):
    return (
        matrix[0][0] * i + matrix[0][1] * j + matrix[0][2] * k + matrix[0][3],
        matrix[1][0] * i + matrix[1][1] * j + matrix[1][2] * k + matrix[1][3],
        matrix[2][0] * i + matrix[2][1] * j + matrix[2][2] * k + matrix[2][3],
    )


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
    values = [float(value) for value in numbers[:9]]
    return [values[0:3], values[3:6], values[6:9]]


def _node_attribute_values(node) -> Iterable[Tuple[str, str]]:
    names = vtk.vtkStringArray()
    if not hasattr(node, "GetAttributeNames"):
        return
    node.GetAttributeNames(names)
    for index in range(names.GetNumberOfValues()):
        key = names.GetValue(index)
        value = node.GetAttribute(key)
        if value:
            yield key, value


def measurement_frame_from_node(volume_node) -> Optional[List[List[float]]]:
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
                    return matrix
        for key, value in _node_attribute_values(node):
            if "measurement" in key.lower() and "frame" in key.lower():
                matrix = _parse_matrix_from_text(value)
                if matrix is not None:
                    return matrix
    return None


def _axis_ranges(dims_ijk: Tuple[int, int, int], settings: GlyphSettings):
    i_min, i_max = 0, dims_ijk[0] - 1
    j_min, j_max = 0, dims_ijk[1] - 1
    k_min, k_max = 0, dims_ijk[2] - 1

    if settings.roi_ijk is not None:
        ri_min, ri_max, rj_min, rj_max, rk_min, rk_max = settings.roi_ijk
        i_min, i_max = max(i_min, ri_min), min(i_max, ri_max)
        j_min, j_max = max(j_min, rj_min), min(j_max, rj_max)
        k_min, k_max = max(k_min, rk_min), min(k_max, rk_max)

    if i_min > i_max or j_min > j_max or k_min > k_max:
        raise ValueError("The requested IJK ROI does not overlap the input volume.")

    if settings.slice_axis != "all":
        if settings.slice_axis == "i":
            index = settings.slice_index if settings.slice_index is not None else (i_min + i_max) // 2
            i_min = i_max = min(max(index, i_min), i_max)
        elif settings.slice_axis == "j":
            index = settings.slice_index if settings.slice_index is not None else (j_min + j_max) // 2
            j_min = j_max = min(max(index, j_min), j_max)
        else:
            index = settings.slice_index if settings.slice_index is not None else (k_min + k_max) // 2
            k_min = k_max = min(max(index, k_min), k_max)

    stride = max(settings.stride, 1)
    return (
        range(i_min, i_max + 1, stride),
        range(j_min, j_max + 1, stride),
        range(k_min, k_max + 1, stride),
    )


def collect_glyph_records(volume_node, mask_node, settings: GlyphSettings):
    warnings: List[str] = []
    validate_input_volume(volume_node)
    array = volume_array_components_last(volume_node)
    dims_ijk = volume_node.GetImageData().GetDimensions()
    expected_kji = (dims_ijk[2], dims_ijk[1], dims_ijk[0])
    if array.shape[:3] != expected_kji:
        raise ValueError(f"Array shape {array.shape[:3]} does not match image dimensions {dims_ijk}.")

    mask_array = mask_array_kji(mask_node) if mask_node is not None else None
    if mask_array is not None and mask_array.shape != expected_kji:
        raise ValueError(f"Mask shape {mask_array.shape} does not match input shape {expected_kji}.")

    measurement_frame = None
    if settings.measurement_frame_mode != "ignore":
        measurement_frame = measurement_frame_from_node(volume_node)
        if measurement_frame is None:
            message = (
                "No NRRD measurement-frame metadata found. Stored directions were assumed to be "
                "already in RAS/physical coordinates."
            )
            if settings.measurement_frame_mode == "require":
                raise ValueError(message)
            warnings.append(message)

    ijk_to_ras = _ijk_to_ras_matrix(volume_node)
    i_range, j_range, k_range = _axis_ranges(dims_ijk, settings)

    records: List[GlyphRecord] = []
    for k in k_range:
        for j in j_range:
            for i in i_range:
                if mask_array is not None and mask_array[k, j, i] <= 0:
                    continue

                center = _point_ijk_to_ras(ijk_to_ras, i, j, k)
                for axis_id, component_start in enumerate(AXIS_COMPONENT_STARTS[: settings.max_axes]):
                    dx, dy, dz, score = array[k, j, i, component_start : component_start + 4]
                    score = float(score)
                    if score <= settings.threshold or not math.isfinite(score):
                        continue

                    direction = _transform_direction((float(dx), float(dy), float(dz)), measurement_frame)
                    if direction is None:
                        continue

                    records.append(GlyphRecord(i, j, k, center, direction, score, axis_id))
                    if len(records) >= settings.max_glyphs:
                        warnings.append(
                            f"Reached the maximum of {settings.max_glyphs} glyphs. "
                            "Increase the limit or use stricter filtering."
                        )
                        return records, warnings
    return records, warnings


def _record_endpoints(record: GlyphRecord, glyph_length: float):
    half = 0.5 * glyph_length
    cx, cy, cz = record.center_ras
    dx, dy, dz = record.direction_ras
    return (
        (cx - dx * half, cy - dy * half, cz - dz * half),
        (cx + dx * half, cy + dy * half, cz + dz * half),
    )


def records_to_polydata(records: Sequence[GlyphRecord], glyph_length: float, tube_radius: float, tube_sides: int):
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
        point0, point1 = _record_endpoints(record, glyph_length)
        point_id0 = points.InsertNextPoint(point0)
        point_id1 = points.InsertNextPoint(point1)
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, point_id0)
        line.GetPointIds().SetId(1, point_id1)
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
        tubed_polydata = vtk.vtkPolyData()
        tubed_polydata.DeepCopy(tube.GetOutput())
        return tubed_polydata
    return polydata


def save_polydata(polydata, export_vtp: str = "", export_vtk: str = "") -> None:
    if export_vtp:
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(export_vtp)
        writer.SetInputData(polydata)
        if writer.Write() != 1:
            raise RuntimeError(f"Failed to write VTP: {export_vtp}")
    if export_vtk:
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(export_vtk)
        writer.SetInputData(polydata)
        if writer.Write() != 1:
            raise RuntimeError(f"Failed to write VTK: {export_vtk}")


def save_csv(records: Sequence[GlyphRecord], export_csv: str) -> None:
    if not export_csv:
        return
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


def remove_generated_models() -> None:
    for model_node in list(slicer.util.getNodesByClass("vtkMRMLModelNode")):
        if model_node.GetAttribute(GENERATED_MODEL_ATTRIBUTE) == "1":
            display_node = model_node.GetDisplayNode()
            if display_node is not None:
                slicer.mrmlScene.RemoveNode(display_node)
            slicer.mrmlScene.RemoveNode(model_node)


def _add_model_node(polydata, name: str, color: Tuple[float, float, float], transform_node_id: Optional[str]):
    model_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", name)
    model_node.SetAttribute(GENERATED_MODEL_ATTRIBUTE, "1")
    model_node.SetAndObservePolyData(polydata)
    if transform_node_id:
        model_node.SetAndObserveTransformNodeID(transform_node_id)
    model_node.CreateDefaultDisplayNodes()
    display_node = model_node.GetDisplayNode()
    display_node.SetColor(color)
    display_node.SetVisibility(True)
    return model_node


def display_records(
    records: Sequence[GlyphRecord],
    glyph_length: float,
    tube_radius: float,
    tube_sides: int,
    model_name: str,
    color_by: str,
    transform_node_id: Optional[str],
):
    model_nodes = []
    if not records:
        return model_nodes

    if color_by == "axis":
        axis_colors = ((0.95, 0.15, 0.15), (0.10, 0.70, 0.25), (0.20, 0.35, 1.00))
        for axis_id, color in enumerate(axis_colors):
            axis_records = [record for record in records if record.axis_id == axis_id]
            if axis_records:
                polydata = records_to_polydata(axis_records, glyph_length, tube_radius, tube_sides)
                model_nodes.append(
                    _add_model_node(polydata, f"{model_name} axis {axis_id + 1}", color, transform_node_id)
                )
    else:
        polydata = records_to_polydata(records, glyph_length, tube_radius, tube_sides)
        model_node = _add_model_node(polydata, model_name, (1.0, 0.85, 0.15), transform_node_id)
        display_node = model_node.GetDisplayNode()
        scores = [record.score for record in records]
        display_node.SetScalarVisibility(True)
        display_node.SetActiveScalarName("score")
        display_node.SetScalarRange(min(scores), max(scores))
        model_nodes.append(model_node)
    return model_nodes


def create_visualization(
    volume_node,
    mask_node,
    settings: GlyphSettings,
    model_name: str,
    export_vtp: str = "",
    export_vtk: str = "",
    export_csv: str = "",
) -> VisualizationResult:
    validate_input_volume(volume_node)
    spacing = volume_node.GetSpacing()
    glyph_length = settings.glyph_length if settings.glyph_length is not None else 0.8 * min(spacing)
    tube_radius = settings.tube_radius if settings.tube_radius is not None else 0.025 * glyph_length

    records, warnings = collect_glyph_records(volume_node, mask_node, settings)
    combined_polydata = records_to_polydata(records, glyph_length, tube_radius, settings.tube_sides)
    save_polydata(combined_polydata, export_vtp, export_vtk)
    save_csv(records, export_csv)

    if settings.clear_previous_models:
        remove_generated_models()
    transform_node_id = volume_node.GetTransformNodeID()
    model_nodes = display_records(
        records,
        glyph_length,
        tube_radius,
        settings.tube_sides,
        model_name,
        settings.color_by,
        transform_node_id,
    )
    return VisualizationResult(len(records), model_nodes, glyph_length, tube_radius, warnings)
