#ifndef SOLVERS_H
#define SOLVERS_H

#include "rdgls_types.h"

//precisionType - pT, MatrixType - MT, VectorType - VT
template <class pT, class MT, class VT>
class SOLVERS
{
public:
	MT *A; //rigelets basis
	MT *s; //input voxel(s)
	pT lmd; //labmda parameter for FISTA training

	SOLVERS();
	~SOLVERS();

	SOLVERS(MT& ridgelets, MT& voxels);
	SOLVERS(MT& ridgelets, MT& voxels, pT lambda);
	SOLVERS(MT& ridgelets, VT& voxel);
	SOLVERS(MT& ridgelets, VT& voxel, pT lambda);

	void FISTA(MT & x, int N_splits);
	void FISTA(VT & x);
};

#include "../src/SOLVERS.cpp"

#endif // !SOLVERS
