[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5591085.svg)](https://doi.org/10.5281/zenodo.5591085)

# Overview

Package to compute spherical ridgelets for diffusion MRI (dMRI).

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# Build from source

The recommended build is the standalone command-line application. If `Eigen3_DIR` and `ITK_DIR` are empty or not provided, CMake downloads and builds pinned compatible versions of Eigen and ITK automatically.

## Requirements

- CMake 3.22 or newer
- Git
- A C++17 compiler
- Ninja or Make
- OpenMP support is recommended for better performance

The first build can take a while because ITK is compiled from source.

## Clone

```sh
git clone https://github.com/rmukh/spherical_ridgelets.git
cd spherical_ridgelets
```

## Windows

The recommended Windows toolchain is MSYS2 UCRT64.

1. Install MSYS2 from <https://www.msys2.org/>.
2. Open the **MSYS2 UCRT64** shell.
3. Install the build tools:

```sh
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja git
```

Then build from the repository root:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1
cmake --build build
```

The executable is created at:

```text
build/SphericalRidgeletsStandalone-build/sphridg.exe
```

MinGW executables need the MSYS2 runtime DLLs available at run time. The build copies the known required DLLs next to `sphridg.exe`. Running from an MSYS2 UCRT64 shell, or adding `C:\msys64\ucrt64\bin` to `PATH`, also works.

If you build from PowerShell or Command Prompt instead of the UCRT64 shell, pass the toolchain explicitly:

```powershell
C:\msys64\ucrt64\bin\cmake.exe -S . -B build -G Ninja `
  -DCMAKE_MAKE_PROGRAM=C:/msys64/ucrt64/bin/ninja.exe `
  -DCMAKE_C_COMPILER=C:/msys64/ucrt64/bin/gcc.exe `
  -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/g++.exe `
  -DJUST_BUILD=1
C:\msys64\ucrt64\bin\cmake.exe --build build
```

## Linux

Install a compiler, CMake, Ninja, and Git. For example:

```sh
# Ubuntu/Debian
sudo apt update
sudo apt install build-essential cmake ninja-build git

# Fedora
sudo dnf install gcc gcc-c++ cmake ninja-build git

# Arch Linux
sudo pacman -S base-devel cmake ninja git
```

If your distribution provides CMake older than 3.22, install a newer CMake from your package manager backports or from <https://cmake.org/download/>.

Build from the repository root:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1
cmake --build build
```

The executable is created at:

```text
build/SphericalRidgeletsStandalone-build/sphridg
```

## macOS

Install the command-line tools and Homebrew packages:

```sh
xcode-select --install
brew install cmake ninja git
```

Build from the repository root:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1
cmake --build build
```

The executable is created at:

```text
build/SphericalRidgeletsStandalone-build/sphridg
```

Apple Clang may not provide OpenMP out of the box. The project can build without OpenMP, but it will be slower. To build with OpenMP support, install LLVM and libomp, then configure with Homebrew LLVM:

```sh
brew install llvm libomp
cmake -S . -B build -G Ninja -DJUST_BUILD=1 \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build
```

## Build options

CMake defaults to `Release` mode. To build in `Debug` mode:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1 -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

To build with float (single precision) instead of double precision:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1 -DUSE_FLOAT=1
cmake --build build
```

To force CMake to ignore any existing `Eigen3_DIR` or `ITK_DIR` values and download the pinned dependencies:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1 -DSPH_DOWNLOAD_DEPS=ON
cmake --build build
```

## Using prebuilt Eigen and ITK

This package mainly depends on Eigen and ITK. If you already have compatible versions installed, pass their CMake package directories:

```sh
cmake -S . -B build -G Ninja -DJUST_BUILD=1 \
  -DEigen3_DIR=/path/to/eigen/share \
  -DITK_DIR=/path/to/ITK-build
cmake --build build
```

`Eigen3_DIR` must contain `Eigen3Config.cmake` or `eigen3-config.cmake`. `ITK_DIR` must contain `ITKConfig.cmake`. ITK 5 or newer is expected; ITK 4.x is rejected or ignored by this build.

Building with `-DJUST_BUILD=1` is essential for a standalone application. Without it, CMake includes the package as a library.

# Basic usage

Mandatory input argument:

- `-i` [dMRI file name]

Optional input arguments:

- `-m` [Mask file]
- `-lvl` [Icosahedron tessellation order, 4 by default]
- `-nspl` [The number of ridgelets coefficient splits for parallel computing, computed automatically by default based on your computer configuration]
- `-mth` [Find maxima ODF threshold, 0.7 by default]
- `-lmd` [Lambda parameter for FISTA solver, 0.01 by default]
- `-sj` [Predefined integer J, which defines the highest level of detectable signal details for the spherical ridgelets, 2 by default]
- `-srho` [Scaling parameter of the spherical ridgelets, 3.125 by default]
- `-rth` [Antipodal non-maximum suppression angle for direct ridgelet maxima, 20 degrees by default]
- `-nth` [The number of threads to use for computations. Otherwise, all available CPU resources will be utilized]
- `-ext_grads` [The external gradients file of (# of directions, 3) shape]
- `-fi` [The number of FISTA iterations, 2000 by default]
- `-ft` [The convergence tolerance of FISTA, 0.001 by default]

Output arguments:

- `-ridg` [Ridgelets coefficients file name]
- `-sr` [Signal reconstruction]
- `-ext_sr` [Signal reconstruction using an external gradients table (`-ext_grads` must be specified)]
- `-odf` [ODF values file name]
- `-omd` [ODF maxima directions and values file name]
- `-omd_r` [Direct ridgelet-coefficient maxima directions and values file name]
- `-A` [A basis file name]
- `-sw` Prints the normalized ridgelet scale weights and exits unless another output is requested
- `-test_omd_r` Runs a synthetic single-coefficient direct ridgelet maxima self-check and exits unless another output is requested
- `-c` Enables compression of output NRRD files, disabled by default

You **must** provide at least one input dMRI file and one output file to run the program, except for `-sw` and `-test_omd_r`, which can be used by themselves.

For example:

```sh
./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd
```

# Outputs

`-ridg` gives a 4D file with the same spatial size as an input.

`-sr` provides a reconstructed signal **without** b0 volumes with the same spatial size as an input.

`-ext_sr` generates a reconstructed signal **without** b0 volumes at the diffusion-encoding directions shipped with the external gradient table.

`-A` outputs a spherical ridgelets basis.

**IMPORTANT!** The output of reconstructed images is always saved with the first dimension representing diffusion-encoding directions.

# Notes on ODF and its directions

The output file for the ODF maximum directions (`-omd`) has the shape of the input dMRI file. Each voxel contains ODF directions and ODF values organized as (x, y, z, ODF value) for each direction. The maximum number of directions is currently fixed at 6 (3 directions, each with an antipode).

The direct ridgelet maxima output (`-omd_r`) uses the same 6-entry `(x, y, z, value)` layout, but estimates directions directly in coefficient space instead of evaluating the full ODF grid. The score for each ridgelet atom is its positive coefficient part weighted by the precomputed aligned scale response. Antipodal non-maximum suppression uses `abs(dot(a, b))` because `v` and `-v` represent the same fiber axis.

For a quick developer check of the scale design:

```sh
./sphridg -sw
```

The default `J = 2` and `srho = 3.125` should print 3 finite scale weights. A synthetic single-coefficient check is available with:

```sh
./sphridg -test_omd_r
```

It sets one finest-scale ridgelet coefficient to 1 and checks that the direct estimator returns that atom's native orientation up to antipodal sign. To compare the direct estimator with the ODF-grid path on a small image, run both `-omd` and `-omd_r` and compare axes up to antipodal sign; exact agreement is not expected because `-omd` peak-picks a sampled ODF while `-omd_r` works directly on native ridgelet orientation grids.

# Important notes

Pre-normalized (by b0) images with no b0 volumes are supported.

If you save NRRD output with an **external** diffusion-encoding directions file, those directions are saved in the metadata, and the original gradients are overridden.

Currently, only the NRRD file format (`.nrrd`, `.nhdr`) is supported.

The input diffusion MRI image is expected to have the shape (size x, size y, size z, # of gradient directions), while the mask file is expected to have the shape (size x, size y, size z, 1). The external gradient file, if used, should not contain comments and should start from the first line as a (# directions, 3) ASCII file.

Saving ODF values **might fail** if you do not have enough RAM.

Additional information is available here: <https://rinatm.com/spherical-ridgelets-for-high-angular-resolution-diffusion-imaging-hardi-implementation/>

TODO: Advanced text cleaning and gradient table detection procedures.

# Advanced users

## Speed

OpenMP is enabled by default when available and provides significant acceleration. If you do not have OpenMP support, we recommend using a compiler with OpenMP support. The split parameter (`-nspl`) is computed to enable a high level of parallelization, although tests have primarily been made on Intel CPUs. If the default value is not optimal for your case, you are encouraged to experiment with it. You can also increase `-nspl` to reduce RAM usage.

Advanced CPU features, such as SSE and AVX, are enabled by default in Release builds when supported by the compiler.

# Bugs and features

Please open a new issue or send a pull request if you found a bug or want to propose a feature. Describe the modification or improvement so it can be reviewed efficiently.

# Citation

**BibTeX**

```bibtex
@misc{sphridg_software,
  title={Software for computation of spherical ridgelets for diffusion MRI},
  DOI={10.5281/zenodo.5591084},
  abstractNote={C++ Package to compute spherical ridgelets.},
  publisher={Zenodo},
  author={Rinat Mukhometzianov and Oleg Michailovich and Yogesh Rathi},
  year={2021},
  month={Oct}
}
```
