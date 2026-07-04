"""Native 3D Slicer GUI for direct spherical-ridgelet direction glyphs."""

import logging

import qt
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleWidget,
)

from RidgeletDirectionsLib import GlyphSettings, create_visualization, validate_input_volume


class RidgeletDirections(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Ridgelet Directions"
        self.parent.categories = ["Diffusion"]
        self.parent.dependencies = []
        self.parent.contributors = ["Rinat Mukhometzianov"]
        self.parent.helpText = (
            "Visualize and export direct spherical-ridgelet direction estimates from "
            "24-component -omd_r NRRD volumes."
        )
        self.parent.acknowledgementText = "Built for the spherical_ridgelets project."


class RidgeletDirectionsWidget(ScriptedLoadableModuleWidget):
    SLICE_AXIS_VALUES = {"I": "i", "J": "j", "K": "k", "All volume": "all"}

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = RidgeletDirectionsLogic()
        self._buildInputsSection()
        self._buildFilteringSection()
        self._buildAppearanceSection()
        self._buildExportSection()
        self._buildApplySection()
        self.layout.addStretch(1)
        self._connectSignals()
        self._updateApplyButton()

    def cleanup(self):
        pass

    def _newNodeSelector(self, node_types):
        selector = slicer.qMRMLNodeComboBox()
        selector.nodeTypes = node_types
        selector.noneEnabled = True
        selector.addEnabled = False
        selector.removeEnabled = False
        selector.renameEnabled = False
        selector.setMRMLScene(slicer.mrmlScene)
        return selector

    def _newLoadButton(self, text, callback):
        button = qt.QPushButton(text)
        button.toolTip = "Load a volume from disk and select it."
        button.connect("clicked(bool)", callback)
        return button

    def _nodeRow(self, selector, button):
        widget = qt.QWidget()
        layout = qt.QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(selector, 1)
        layout.addWidget(button)
        return widget

    def _buildInputsSection(self):
        section = qt.QGroupBox("Inputs")
        form = qt.QFormLayout(section)

        self.omdSelector = self._newNodeSelector(
            ["vtkMRMLVectorVolumeNode", "vtkMRMLScalarVolumeNode"]
        )
        self.referenceSelector = self._newNodeSelector(
            ["vtkMRMLScalarVolumeNode", "vtkMRMLVectorVolumeNode"]
        )
        self.maskSelector = self._newNodeSelector(
            ["vtkMRMLLabelMapVolumeNode", "vtkMRMLScalarVolumeNode"]
        )

        form.addRow(
            "Direct directions:",
            self._nodeRow(self.omdSelector, self._newLoadButton("Load -omd_r NRRD...", self.onLoadOmd)),
        )
        form.addRow(
            "Reference volume:",
            self._nodeRow(
                self.referenceSelector,
                self._newLoadButton("Load reference...", self.onLoadReference),
            ),
        )
        form.addRow(
            "Mask volume:",
            self._nodeRow(self.maskSelector, self._newLoadButton("Load mask...", self.onLoadMask)),
        )
        self.layout.addWidget(section)

    def _buildFilteringSection(self):
        section = qt.QGroupBox("Voxel filtering")
        form = qt.QFormLayout(section)

        self.sliceAxisCombo = qt.QComboBox()
        self.sliceAxisCombo.addItems(["K", "J", "I", "All volume"])
        self.sliceIndexSpin = qt.QSpinBox()
        self.sliceIndexSpin.setRange(0, 0)
        slice_row = qt.QWidget()
        slice_layout = qt.QHBoxLayout(slice_row)
        slice_layout.setContentsMargins(0, 0, 0, 0)
        slice_layout.addWidget(self.sliceAxisCombo)
        slice_layout.addWidget(self.sliceIndexSpin)
        form.addRow("Slice axis / index:", slice_row)

        self.strideSpin = qt.QSpinBox()
        self.strideSpin.setRange(1, 100)
        self.strideSpin.setValue(2)
        form.addRow("Stride:", self.strideSpin)

        self.thresholdSpin = qt.QDoubleSpinBox()
        self.thresholdSpin.setDecimals(6)
        self.thresholdSpin.setRange(-1.0e12, 1.0e12)
        self.thresholdSpin.setValue(0.0)
        form.addRow("Minimum score:", self.thresholdSpin)

        self.maxAxesSpin = qt.QSpinBox()
        self.maxAxesSpin.setRange(1, 3)
        self.maxAxesSpin.setValue(3)
        form.addRow("Maximum axes per voxel:", self.maxAxesSpin)

        self.maxGlyphsSpin = qt.QSpinBox()
        self.maxGlyphsSpin.setRange(1, 10000000)
        self.maxGlyphsSpin.setValue(100000)
        self.maxGlyphsSpin.setSingleStep(10000)
        form.addRow("Maximum glyphs:", self.maxGlyphsSpin)

        self.useRoiCheck = qt.QCheckBox("Restrict to inclusive IJK bounds")
        form.addRow("IJK ROI:", self.useRoiCheck)
        self.roiSpins = []
        roi_grid = qt.QGridLayout()
        for column, axis in enumerate(("I", "J", "K")):
            roi_grid.addWidget(qt.QLabel(axis + " min"), 0, column * 2)
            minimum = qt.QSpinBox()
            minimum.setRange(0, 0)
            roi_grid.addWidget(minimum, 0, column * 2 + 1)
            roi_grid.addWidget(qt.QLabel(axis + " max"), 1, column * 2)
            maximum = qt.QSpinBox()
            maximum.setRange(0, 0)
            roi_grid.addWidget(maximum, 1, column * 2 + 1)
            self.roiSpins.extend((minimum, maximum))
        self.roiWidget = qt.QWidget()
        self.roiWidget.setLayout(roi_grid)
        self.roiWidget.enabled = False
        form.addRow("", self.roiWidget)

        self.measurementFrameCombo = qt.QComboBox()
        self.measurementFrameCombo.addItems(["Auto", "Ignore", "Require"])
        form.addRow("Measurement frame:", self.measurementFrameCombo)

        self.layout.addWidget(section)

    def _buildAppearanceSection(self):
        section = qt.QGroupBox("Glyph appearance")
        form = qt.QFormLayout(section)

        self.autoLengthCheck = qt.QCheckBox("Automatic from voxel spacing")
        self.autoLengthCheck.checked = True
        self.glyphLengthSpin = qt.QDoubleSpinBox()
        self.glyphLengthSpin.setDecimals(4)
        self.glyphLengthSpin.setRange(0.0001, 100000.0)
        self.glyphLengthSpin.setValue(2.0)
        self.glyphLengthSpin.enabled = False
        length_row = qt.QWidget()
        length_layout = qt.QHBoxLayout(length_row)
        length_layout.setContentsMargins(0, 0, 0, 0)
        length_layout.addWidget(self.autoLengthCheck)
        length_layout.addWidget(self.glyphLengthSpin)
        form.addRow("Glyph length:", length_row)

        self.autoRadiusCheck = qt.QCheckBox("Automatic from glyph length")
        self.autoRadiusCheck.checked = True
        self.tubeRadiusSpin = qt.QDoubleSpinBox()
        self.tubeRadiusSpin.setDecimals(5)
        self.tubeRadiusSpin.setRange(0.0, 10000.0)
        self.tubeRadiusSpin.setValue(0.05)
        self.tubeRadiusSpin.enabled = False
        radius_row = qt.QWidget()
        radius_layout = qt.QHBoxLayout(radius_row)
        radius_layout.setContentsMargins(0, 0, 0, 0)
        radius_layout.addWidget(self.autoRadiusCheck)
        radius_layout.addWidget(self.tubeRadiusSpin)
        form.addRow("Tube radius:", radius_row)

        self.tubeSidesSpin = qt.QSpinBox()
        self.tubeSidesSpin.setRange(3, 32)
        self.tubeSidesSpin.setValue(8)
        form.addRow("Tube sides:", self.tubeSidesSpin)

        self.colorByCombo = qt.QComboBox()
        self.colorByCombo.addItems(["Axis rank", "Score"])
        form.addRow("Color by:", self.colorByCombo)

        self.modelNameEdit = qt.QLineEdit("Direct ridgelet axes")
        form.addRow("Model name:", self.modelNameEdit)

        self.clearPreviousCheck = qt.QCheckBox("Remove models previously created by this module")
        self.clearPreviousCheck.checked = True
        form.addRow("", self.clearPreviousCheck)

        self.layout.addWidget(section)

    def _newExportRow(self, extension, callback):
        edit = qt.QLineEdit()
        edit.placeholderText = f"Optional {extension.upper()} output path"
        button = qt.QPushButton("Browse...")
        button.connect("clicked(bool)", callback)
        row = qt.QWidget()
        layout = qt.QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit, 1)
        layout.addWidget(button)
        return row, edit

    def _buildExportSection(self):
        section = qt.QGroupBox("Optional conversion exports")
        form = qt.QFormLayout(section)
        vtp_row, self.exportVtpEdit = self._newExportRow("vtp", self.onBrowseVtp)
        vtk_row, self.exportVtkEdit = self._newExportRow("vtk", self.onBrowseVtk)
        csv_row, self.exportCsvEdit = self._newExportRow("csv", self.onBrowseCsv)
        form.addRow("VTK XML PolyData:", vtp_row)
        form.addRow("Legacy VTK PolyData:", vtk_row)
        form.addRow("CSV:", csv_row)
        self.layout.addWidget(section)

    def _buildApplySection(self):
        self.applyButton = qt.QPushButton("Create visualization")
        self.applyButton.toolTip = "Generate Slicer model nodes and optional export files."
        self.applyButton.setMinimumHeight(36)
        self.statusLabel = qt.QLabel("Select a direct directions volume.")
        self.statusLabel.wordWrap = True
        self.layout.addWidget(self.applyButton)
        self.layout.addWidget(self.statusLabel)

    def _connectSignals(self):
        self.omdSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onInputChanged)
        self.sliceAxisCombo.connect("currentIndexChanged(int)", self.onSliceAxisChanged)
        self.useRoiCheck.connect("toggled(bool)", self.roiWidget.setEnabled)
        self.autoLengthCheck.connect("toggled(bool)", lambda checked: self.glyphLengthSpin.setDisabled(checked))
        self.autoRadiusCheck.connect("toggled(bool)", lambda checked: self.tubeRadiusSpin.setDisabled(checked))
        self.applyButton.connect("clicked(bool)", self.onApply)

    def _loadVolumeFile(self, title, vector=False):
        path = qt.QFileDialog.getOpenFileName(
            slicer.util.mainWindow(),
            title,
            "",
            "NRRD volumes (*.nrrd *.nhdr);;All files (*)",
        )
        if not path:
            return None
        if vector:
            try:
                loaded = slicer.util.loadNodeFromFile(path, "VectorVolumeFile", {}, True)
                if isinstance(loaded, tuple):
                    success, node = loaded
                    if success:
                        return node
                elif loaded is not None:
                    return loaded
            except Exception as error:
                logging.warning("Vector volume load failed, using generic loader: %s", error)
        loaded = slicer.util.loadVolume(path, returnNode=True)
        if isinstance(loaded, tuple):
            success, node = loaded
            if success:
                return node
        elif loaded is not None:
            return loaded
        raise RuntimeError(f"Could not load volume: {path}")

    def onLoadOmd(self):
        node = self._loadVolumeFile("Load direct ridgelet directions", vector=True)
        if node is not None:
            self.omdSelector.setCurrentNode(node)

    def onLoadReference(self):
        node = self._loadVolumeFile("Load reference volume")
        if node is not None:
            self.referenceSelector.setCurrentNode(node)

    def onLoadMask(self):
        node = self._loadVolumeFile("Load mask volume")
        if node is not None:
            self.maskSelector.setCurrentNode(node)

    def _browseOutput(self, title, file_filter, edit):
        path = qt.QFileDialog.getSaveFileName(slicer.util.mainWindow(), title, "", file_filter)
        if path:
            edit.text = path

    def onBrowseVtp(self):
        self._browseOutput("Export VTP", "VTK XML PolyData (*.vtp)", self.exportVtpEdit)

    def onBrowseVtk(self):
        self._browseOutput("Export VTK", "Legacy VTK PolyData (*.vtk)", self.exportVtkEdit)

    def onBrowseCsv(self):
        self._browseOutput("Export CSV", "CSV files (*.csv)", self.exportCsvEdit)

    def onInputChanged(self, node=None):
        self._updateVolumeRanges()
        self._updateApplyButton()

    def onSliceAxisChanged(self, index=None):
        self._updateSliceIndexRange()

    def _updateApplyButton(self):
        valid = self.omdSelector.currentNode() is not None
        self.applyButton.enabled = valid
        if not valid:
            self.statusLabel.text = "Select a direct directions volume."

    def _updateVolumeRanges(self):
        node = self.omdSelector.currentNode()
        if node is None or node.GetImageData() is None:
            return
        dims = node.GetImageData().GetDimensions()
        for axis, (minimum_spin, maximum_spin) in enumerate(
            ((self.roiSpins[0], self.roiSpins[1]),
             (self.roiSpins[2], self.roiSpins[3]),
             (self.roiSpins[4], self.roiSpins[5]))
        ):
            maximum = max(dims[axis] - 1, 0)
            minimum_spin.setRange(0, maximum)
            maximum_spin.setRange(0, maximum)
            minimum_spin.setValue(0)
            maximum_spin.setValue(maximum)
        self._updateSliceIndexRange()

    def _updateSliceIndexRange(self):
        node = self.omdSelector.currentNode()
        axis_text = self.sliceAxisCombo.currentText
        if node is None or node.GetImageData() is None or axis_text == "All volume":
            self.sliceIndexSpin.enabled = False
            return
        self.sliceIndexSpin.enabled = True
        dims = node.GetImageData().GetDimensions()
        axis_index = {"I": 0, "J": 1, "K": 2}[axis_text]
        maximum = max(dims[axis_index] - 1, 0)
        self.sliceIndexSpin.setRange(0, maximum)
        self.sliceIndexSpin.setValue(maximum // 2)

    def _settingsFromGui(self):
        roi = None
        if self.useRoiCheck.checked:
            roi = tuple(spin.value for spin in self.roiSpins)
        return GlyphSettings(
            slice_axis=self.SLICE_AXIS_VALUES[self.sliceAxisCombo.currentText],
            slice_index=None if self.sliceAxisCombo.currentText == "All volume" else self.sliceIndexSpin.value,
            roi_ijk=roi,
            stride=self.strideSpin.value,
            threshold=self.thresholdSpin.value,
            max_axes=self.maxAxesSpin.value,
            max_glyphs=self.maxGlyphsSpin.value,
            glyph_length=None if self.autoLengthCheck.checked else self.glyphLengthSpin.value,
            tube_radius=None if self.autoRadiusCheck.checked else self.tubeRadiusSpin.value,
            tube_sides=self.tubeSidesSpin.value,
            color_by="axis" if self.colorByCombo.currentText == "Axis rank" else "score",
            measurement_frame_mode=self.measurementFrameCombo.currentText.lower(),
            clear_previous_models=self.clearPreviousCheck.checked,
        )

    def onApply(self):
        with slicer.util.tryWithErrorDisplay("Failed to create ridgelet direction visualization.", waitCursor=True):
            volume_node = self.omdSelector.currentNode()
            validate_input_volume(volume_node)
            reference_node = self.referenceSelector.currentNode()
            if reference_node is not None:
                slicer.util.setSliceViewerLayers(background=reference_node)

            result = self.logic.process(
                volume_node,
                self.maskSelector.currentNode(),
                self._settingsFromGui(),
                self.modelNameEdit.text.strip() or "Direct ridgelet axes",
                self.exportVtpEdit.text.strip(),
                self.exportVtkEdit.text.strip(),
                self.exportCsvEdit.text.strip(),
            )
            self.statusLabel.text = (
                f"Created {len(result.model_nodes)} model node(s) from {result.glyph_count} glyphs. "
                f"Length={result.glyph_length:.4g}, radius={result.tube_radius:.4g}."
            )
            if result.warnings:
                slicer.util.warningDisplay("\n\n".join(result.warnings), windowTitle="Ridgelet Directions")


class RidgeletDirectionsLogic(ScriptedLoadableModuleLogic):
    def process(
        self,
        volume_node,
        mask_node,
        settings,
        model_name,
        export_vtp="",
        export_vtk="",
        export_csv="",
    ):
        return create_visualization(
            volume_node,
            mask_node,
            settings,
            model_name,
            export_vtp,
            export_vtk,
            export_csv,
        )
