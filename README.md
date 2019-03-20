# Overview
Package to compute spherical ridgelets.

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## Easiest way to build (Linux)
1. Download or clone repository

        git clone https://github.com/research-enthusiast/spherical_ridgelets.git

2. Create empty folder inside downloaded repository.
    
        cd spherical_ridgelets
        mkdir build

3. Enter that directory

        cd build

4. Generate make files

        cmake ..

5. Build package

        make

The final binary file will be in *Ridgelets-build* directory named as *sphridg*

## Basic Usage
* Mandatory input argument: *-i [dMRI file name]*
* Optional input arguments: *-m [mask file]*
* Output argumets: *-ridg [ridgelet file name]* -odf *[ODF values file name]* -omd *[ODF maxima directions and values file name]* *-c* enable compression of output files

**Should** be at least at least input dMRI file and one output file.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Important notes
For now this software supports NRRD file formats only (.nrrd, .nhdr) both for input and output. To build this project you should have **CMake** and **git** installed. The repository contains **Visual Studio 2017** project files for current development purposes. Currently **gcc** compiler adequtely supported. Possibly **clang** works fine, but have not been tested yet.

# Advanced users

## Speed
All cmake files created in a way that during the building cmake will automatically determine if you have OpenMP installed and use it during compilation. This gives significant speedup. Also, it is determines and build package with supported CPU features like SSE, AVX, etc.

## Custom build
This package mainly depends on two libraries: ITK and Eigen. In some cases you may want to use custom versions of these libraries. That's typically happens in the following cases:
* You have them preinstalled in your system and want to save time on compilation process 
* You want to test some specific version for perfomace comparasion
than you can pass path in cmake command.
* For ITK pass -DITK_DIR
* For Eigen pass -DEigen3_DIR

For example:

        cmake -DITK_DIR=/path/to/ITK_build -DEigen3_DIR=/path/to/Eigen ..
