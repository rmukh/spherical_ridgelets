# Overview
Library to compute spherical ridgelets.

Authors: Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

# How to start

## Easiest way to build
1. Download or clone repository (git clone https://github.com/restearch-enthusiast/spherical_ridgelets.git) 
2. Create empty folder inside downloaded repository. *build* in that example
3. Enter that directory *cd build*
4. *cmake ..*
5. *make*

## Basic Usage
./sphridg -i <dMRI file name> and at least one output: -ridg, -odf, -om, -or
  
*Optional input arguments*: -m <mask file>

*Possible output argumet(s)*: -ridg ridgelet_file -odf ODF_values -om ODF_maxima -or orientations -c enable compression

# Notes
For now only NRRD file formats supported (.nrrd, .nhdr)
