# Overview
Package to compute spherical ridgelets.

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## The easiest way to build as a standalone application (Linux)

You can run the script *linux_standalone_build.sh* to perform all steps described below.

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
- *-A* [A basis file]
- *-c* enables compression of output files, false by default

You **must** provide at least input dMRI file and one output file to make any computations possible.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Notes on ODF and its directions
Output file for *-omd* has a shape of input dMRI file and each voxel contains ODF directions and ODF values in that directions organized as (x y z odf_value) for each direction. Now maximum number of directions fixed to 6.

# Important notes
It is very important to build it with flag *-DJUST_BUILD=1*. Otherwise, the CMakeLists.txt will include files necessary to build this package as a library.

For now, this software supports NRRD file formats only (.nrrd, .nhdr) both for input and output. To build this project you should have *CMake* and *git* installed.

Input diffusion MRI image expected to be in the shape of (size_x, size_y, size_z, #_of_gradient directions), while mask file expected to be in the shape of (size_x, size_y, size_z, 1).

The repository contains *Visual Studio 2017* project files for current development purposes. Currently *gcc* compiler adequately supported. Possibly *clang* works fine, but have not been tested yet. This package was tested on *Linux* only! *MacOS* and *Windows* compatibility is not guaranteed for now. When saving **ODF values** you may experience problems with that if you don't have enough **RAM memory** and saving operation might **fail**.

# Advanced users

## Speed
All cmake files created in a way that during the building cmake will automatically determine if you have *OpenMP* installed and use it during compilation. This gives a significant speedup. So, if you don't have it installed, consider its installation. [Google](https://www.google.com/) and [this page](https://www.openmp.org/resources/openmp-compilers-tools/) are excellent sources of information on *OpenMP*. Also, it detects and builds the package with supported CPU features like SSE, AVX, etc.

By default *-nspl* parameter is computed in a way which provides the highest possible level of parallelization for that implementation. It was tested on Intel CPUs. If you feel that the default value is not optimal for your case, feel free to experiment with that parameter. You can increase default value of *-nspl* to reduce RAM consumption.

By default cmake starts building this package and all required libraries in Release mode to achieve fastest perfomance. If you want to build in Debug mode, don't forget to pass -DCMAKE_BUILD_TYPE=Debug

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

Please, open a new issue if you found any bug, error or just want to propose new features. Also, feel free to send pull requests, create new branches and send merge request. Please, provide a description of the modifications or/and improvements you made, otherwise it might take longer for me to review and accept your request.
