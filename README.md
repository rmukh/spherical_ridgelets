# Overview
Package to compute spherical ridgelets.

**Authors**: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## Easiest way to build (Linux)
1. Download or clone repository (*git clone https://github.com/research-enthusiast/spherical_ridgelets.git*) 
2. Create empty folder inside downloaded repository. (*mkdir build*)
3. Enter that directory *cd build*
4. *cmake ..*
5. *make*

The final binary file will be in *Ridgelets-build* directory named as *sphridg*

## Basic Usage
* Mandatory input argument: *-i [dMRI file name]*
* Optional input arguments: *-m [mask file]*
* Output argumets: *-ridg [ridgelet file name]* -odf *[ODF values file name]* -omd *[ODF maxima directions and values file name]* *-c* enable compression of output files

**Should** be at least at least input dMRI file and one output file.

For example:

    ./sphridg -i my_dmri.nrrd -ridg ridgelets_coefficients.nrrd

# Notes
For now this software supports NRRD file formats only (.nrrd, .nhdr) both for input and output. To build this project you should have **CMake** and **git** installed. The repository contains **Visual Studio 2017** project files for current development purposes.
