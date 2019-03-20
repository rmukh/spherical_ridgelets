# Overview
Package to compute spherical ridgelets.

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## The easiest way to build (Linux)
1. Download or clone repository.

        git clone https://github.com/research-enthusiast/spherical_ridgelets.git

2. Create an empty folder inside the downloaded repository.
    
        cd spherical_ridgelets
        mkdir build

3. Enter that directory.

        cd build

4. Generate make files.

        cmake ..

5. Build package.

        make

The final binary file will be in *Ridgelets-build* directory named as *sphridg*

## Basic Usage
* Mandatory input argument: *-i [dMRI file name]*
* Optional input arguments: *-m [mask file]*
* Output arguments: *-ridg [ridgelet file name]* -odf *[ODF values file name]* -omd *[ODF maxima directions and values file name]* *-c* enable compression of output files

**Should** be at least input dMRI file and one output file.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Important notes
For now, this software supports NRRD file formats only (.nrrd, .nhdr) both for input and output. To build this project you should have *CMake* and *git* installed. The repository contains *Visual Studio 2017* project files for current development purposes. Currently *gcc* compiler adequately supported. Possibly *clang* works fine, but have not been tested yet. This package was tested on *Linux* only! *MacOS* and *Windows* compatibility is not guaranteed for now.

# Advanced users

## Speed
All cmake files created in a way that during the building cmake will automatically determine if you have *OpenMP* installed and use it during compilation. This gives a significant speedup. So, if you don't have it installed, consider its installation. [Google](https://www.google.com/) and [this page](https://www.openmp.org/resources/openmp-compilers-tools/) are excellent sources of information on *OpenMP*. Also, it detects and builds the package with supported CPU features like SSE, AVX, etc. 

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
