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

Note: You can also use *make* with multi-threading option in a way:

        make -j#of_threads, e.g. make -j4

The final binary file will be in *build/Spherical_Ridgelets-build* directory and named as *sphridg*

## Basic Usage

Mandatory input argument:
- *-i* [dMRI file name]

Optional input arguments:
- *-m* [Mask file]
- *-lvl* [icosahedron tesselation order, 4 by default]
- *-nspl* [The number of ridgelets coefficients splits to parallel computing, computed automatically by default based on your computer configuration]
- *-mth* [Find maxima ODF threshold, 0.7 by default] 
- *-lmd* [Lambda parameter for FISTA solver, 0.01 by default] 
- *-sj* [Predefined integer J, which defines the highest level of 'detectable' signal details parameter of the spherical ridgelets, 2 by default] 
- *-srho* [Scaling parameter of the spherical ridgelets, 3.125 by default]
- *-nth* [The number of threads to use for computations. Otherwise, all available CPU resources will be utilized]

Output arguments:
- *-ridg* [Ridgelet file name]
- *-sr* [Signal reconstruction]
- *-odf* [ODF values file name] 
- *-omd* [ODF maxima directions and values file name]
- *-A* [A basis file name]
- *-c* enables compression of output files, false by default

You **must** provide at least input dMRI file and one output file to make any computations possible.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Notes on ODF and its directions
Output file for the ODF maximum directions (-omd) has a size of input dMRI file. Each voxel contains ODF directions and ODF values organized as (x y z odf_value) for each direction. Now a maximum number of directions is fixed to 6 (3 directions, each has an antipode).

# Important notes
It is very important to build it with the flag *-DJUST_BUILD=1*. Otherwise, the CMakeLists.txt will include files necessary to build this package as a library.

Currently, NRRD file format (.nrrd, .nhdr) supported only. To build this project you need [CMake](https://cmake.org/) and [git](https://git-scm.com/) in your system. 

Input diffusion MRI image expected to be in the shape of (size_x, size_y, size_z, #_of_gradient directions), while mask file expected to be in the shape of (size_x, size_y, size_z, 1).

The repository contains *Visual Studio 2017* project files for the development purposes, so you can safely delete them. *GCC* compiler adequately supported. Possibly *clang* works fine, but have not been tested yet. This package was tested on *Linux* only! *MacOS* and *Windows* compatibility is not guaranteed. Saving *ODF values* operation might **fail** if you don't have **enough RAM memory**.
# Advanced users

## Speed
All cmake files created in a way that during the building cmake will automatically determine if you have *OpenMP* installed and use it during compilation. This gives a significant speedup. So, if you don't have OpenMP support, we recommend you to use GCC compiler with OpenMP. [Google](https://www.google.com/) and [this page](https://www.openmp.org/resources/openmp-compilers-tools/) are excellent sources of information on *OpenMP*. Also, it detects and builds the package with supported CPU features like SSE, AVX, etc.

The split parameter (*-nspl*) computed in a way to provide the highest possible level of parallelization for that implementation. This parallelization was tested on Intel only. If you feel that the default value is not optimal for your case, feel free to experiment with that parameter. You can also increase the default value of *-nspl* to reduce RAM consumption.

CMake starts building this package and all required libraries in Release mode to achieve the highest performance. If you want to build in Debug mode, don't forget to pass -DCMAKE_BUILD_TYPE=Debug

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

Please, open a new issue or send a pull request if you found any bug/error or want to propose new features. Don't forget to provide a description of the modifications or/and improvements you made, otherwise it might take longer for me to review and accept your request.
