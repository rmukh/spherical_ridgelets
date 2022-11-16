[![Build Status](https://app.travis-ci.com/rmukh/spherical_ridgelets.svg?branch=master)](https://app.travis-ci.com/rmukh/spherical_ridgelets) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5591085.svg)](https://doi.org/10.5281/zenodo.5591085)


# Overview
Package to compute spherical ridgelets.

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## The easiest way is to build a standalone application (Linux)

You can run the script **linux_standalone_build.sh** to perform all steps described below.

1. Download or clone repository.

        git clone https://github.com/rmukh/spherical_ridgelets.git

2. Create an empty folder inside the downloaded repository.
    
        cd spherical_ridgelets
        mkdir build

3. Enter that directory.

        cd build

4. Generate make files.

        cmake .. -DJUST_BUILD=1

5. Build package.

        make

### Notes: 

1. You can also use *make* with multi-threading option in a way:

        make -j#of_threads, e.g. make -j4

2. If you want to compile the package with float (single) precision type used instead of double, please add *-DUSE_FLOAT=1* when calling the CMake tool.

3. Please, find the binary file *sphridg* in the *build/Spherical_Ridgelets-build* directory if building succeeded. 

## Basic Usage

Mandatory input argument:
- *-i* [dMRI file name]

Optional input arguments:
- *-m* [Mask file]
- *-lvl* [Icosahedron tesselation order, 4 by default]
- *-nspl* [The number of ridgelets coefficients splits to parallel computing, computed automatically by default based on your computer configuration]
- *-mth* [Find maxima ODF threshold, 0.7 by default] 
- *-lmd* [Lambda parameter for FISTA solver, 0.01 by default] 
- *-sj* [Predefined integer J, which defines the highest level of 'detectable' signal details parameter of the spherical ridgelets, 2 by default] 
- *-srho* [Scaling parameter of the spherical ridgelets, 3.125 by default]
- *-nth* [The number of threads to use for computations. Otherwise, all available CPU resources will be utilized]
- *-ext_grads* [The external gradients file of (# of directions, 3) shape]
- *-fi* [The number of FISTA iterations, 2000 by default]
- *-ft* [The convergence tolerance of FISTA, 0.001 by default]

Output arguments:
- *-ridg* [Ridgelets coefficients file name]
- *-sr* [Signal reconstruction]
- *-ext_sr* [Signal reconstruction using an external gradients table (*-ext_grads* must be specified)]
- *-odf* [ODF values file name] 
- *-omd* [ODF maxima directions and values file name]
- *-A* [A basis file name]
- *-c* Enables compression of output nrrd files, disabled by default

You **must** provide at least input dMRI file and one output file to run the program.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Outputs
*-ridg* gives a 4D file with the same spatial size as an input, and the last dimension is a vector of representation (ridgelets) coefficients.

*-sr* provides a reconstructed signal **without** b0 volumes with the same spatial size as an input.

*-ext_sr* generates a reconstructed signal **without** b0 volumes at the diffusion-encoding directions shipped with the external gradient table.

*-A* outputs a spherical ridgelets basis.

# Notes on ODF and its directions
Output file for the ODF maximum directions (-omd) has a shape of input dMRI file. Each voxel contains ODF directions and ODF values organized as (x, y, z, odf value) for each direction. Now a maximum number of directions is fixed to 6 (3 directions, each has an antipode).

# Important notes
All **b0** volumes should be **in the beginning** aka first voxels. They couldn't be spread around, located between diffusion-encoded ones, or placed in the end. **Instead**, you can use pre-normalized images with no b0 volumes.

If you are saving NRRD output with an **external** diffusion-encoding directions file, they will be saved in the meta-data information; hence, the original gradients will be overridden.

Building it with the flag *-DJUST_BUILD=1* is essential if you want standalone software. Otherwise, the CMakeLists.txt will include files necessary to make this package in the form of a library.

Currently, the NRRD file format (.nrrd, .nhdr) is supported only. To build this project, you need [CMake](https://cmake.org/) and [git](https://git-scm.com/) installed on your system. 

Input diffusion MRI image expected to be in the shape of (size x, size y, size z, # of gradient directions), while mask file expected to be in the shape of (size x, size y, size z, 1). The external gradient file (if used) should not contain any comments and start from the first line, so just (#directins, 3) ASCII file. The advanced text cleaning and gradient table detection procedures still need to be implemented.

The repository contains *Visual Studio* project files for development purposes, so you can safely delete them. *GCC*, *Clang*, *MSVC* compilers adequately supported. This package was tested on *Linux*, *Windows 10*, *Mac OS*. Please, refer to the Travis CI badge at the top of this manual to check if the current version is compilable. Only *GCC* version 7.x.x is currently adequately supported, so install and use this version (GCC and g++) and specify the system paths if necessary. The usage example for Cmake/make build you can find in linux_standalone_build.sh

Saving ODF values operation **might fail** if you don't have enough RAM.

Addition info you can find here: https://rinatm.com/spherical-ridgelets-for-high-angular-resolution-diffusion-imaging-hardi-implementation/


# Advanced users

## Speed
*OpenMP* enabled by default during compilation to provide you a significant acceleration. So, if you don't have OpenMP support, we recommend you to use GCC compiler with OpenMP. [Google](https://www.google.com/) and [this page](https://www.openmp.org/resources/openmp-compilers-tools/) are excellent sources of information on *OpenMP*. Also, advanced CPU features, e.g., SSE, AVX, etc. enabled by default.

The split parameter (*-nspl*) computed in a way to enable the highest possible level of parallelization, however, tests made on Intel CPU's only. If you feel that the default value is not optimal for your case, you are encouraged to experiment with it. You can also increase the default value of *-nspl* to reduce RAM usage.

CMake starts building this package and all required libraries in Release mode to achieve the highest performance. If you want to build in Debug mode, don't forget to pass *-DCMAKE_BUILD_TYPE=Debug* in CMake

## Custom build
This package mainly depends on two libraries: *ITK* and *Eigen*. In some cases, you may want to use custom versions of these libraries. That's typically happening in the following cases:
* You have them preinstalled in your system and want to save time on the compilation process;
* You want to test some specific version for performance comparasion.

Then you can pass the path in cmake command.
* For *ITK* use *-DITK_DIR*
* For *Eigen* use *-DEigen3_DIR*

For example:

    cmake -DITK_DIR=/path/to/ITK_build -DEigen3_DIR=/path/to/Eigen ..

# Bugs & feautures

Please, open a new issue or send a pull request if you found any bug/error or want to propose new features. Don't forget to describe the modifications or/and improvements you made otherwise it might take a long time to review and accept your request.

# Citation

**BibTeX**  
@misc{sphridg_software, title={Sofwater for computation of spherical ridgelets for diffusion MRI}, DOI={10.5281/zenodo.5591084}, abstractNote={C++ Package to compute spherical ridgelets.}, publisher={Zenodo}, author={Rinat Mukhometzianov and Oleg Michailovich and Yogesh Rathi}, year={2021}, month={Oct} }
