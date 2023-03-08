"""
An examplar file demonstrating how to compute spherical ridgelets coefficients for a 81 voxels.

The npy file located within the 'Demos' folder is obtained from the 'dipy' v.1.6.0 package sample data.
In particular, HARDI150.nii.gz downloaded via fetch_stanford_hardi() function.
The selected voxels are dwi[:, 52, 32, :]

The diffusion encoding voxels (150 directions) are normalized by the mean b0.
"""
import numpy as np
from os.path import join
import Spherical_Ridgelets_py

data = np.load(join('demo81','data.npy'))
bvals = np.load(join('demo81','bvals.npy'))
bvecs = np.load(join('demo81','bvecs.npy'))

print("Data shape", data.shape)

# Fortran ordering
bvecs = np.asfortranarray(bvecs)

# initialize ridgelets basis
sph_ridg = Spherical_Ridgelets_py.SPH_RIDG()
sph_ridg.init()

# Compute ridgelets basis
A = sph_ridg.GetRBasis(bvecs)
print("A shape", A.shape)

# Take a single voxel
signal = np.asfortranarray(data[40,:])
print("signal", signal.shape)

# FISTA solver parameters
n_iterations = 2000
lmd = 0.01
tolerance = 1e-6

# run computations
solver = Spherical_Ridgelets_py.SOLVERS()
C = solver.SolveSingleSignal(A, signal, lmd, n_iterations, tolerance)

print("C shape", C.shape)
print(C)
