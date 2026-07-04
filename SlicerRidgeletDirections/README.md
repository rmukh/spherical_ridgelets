# Slicer Ridgelet Directions

`RidgeletDirections` is a scripted 3D Slicer module for visualizing and
converting the custom 24-component NRRD files produced by spherical ridgelets
`-omd_r`.

Each voxel stores three ranked fiber axes and their antipodes:

```text
dir1, -dir1, dir2, -dir2, dir3, -dir3
```

The module renders the three unique axes as centered tube/stick glyphs.

## Install for development

1. Open 3D Slicer.
2. Open **Edit > Application Settings > Modules**.
3. Add the `SlicerRidgeletDirections/RidgeletDirections` directory to
   **Additional module paths**.
4. Restart Slicer.
5. Open **Modules > Diffusion > Ridgelet Directions**.

The module folder is self-contained and may be moved into its own repository
later.

## Use

1. Load a direct ridgelet directions NRRD with **Load -omd_r NRRD**, or select
   an already-loaded vector volume.
2. Optionally select a reference volume and mask.
3. Choose slice/ROI filtering, stride, score threshold, maximum axes, and
   glyph appearance.
4. Optionally select VTP, legacy VTK, or CSV export paths.
5. Click **Create visualization**.

The default view renders the middle `K` slice at stride 2. Axis colors are:

- primary axis: red
- second axis: green
- third axis: blue

When measurement-frame metadata is unavailable, the module warns and assumes
the stored directions are already in RAS/physical coordinates.

## Build as an extension

Configure this directory against a Slicer build:

```sh
cmake -S . -B build -DSlicer_DIR=<SLICER_BUILD_DIR>
cmake --build build
```
