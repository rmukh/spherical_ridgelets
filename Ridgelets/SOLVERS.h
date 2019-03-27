#ifndef SOLVERS_H
#define SOLVERS_H

#include "UtilMath.h"

class SOLVERS
{
public:
	MatrixType *A; //rigelets basis
	MatrixType *s; //input voxel(s)
	double lmd; //labmda parameter for FISTA training

	SOLVERS();
	~SOLVERS();

	SOLVERS(MatrixType& ridgelets, MatrixType& voxels);
	SOLVERS(MatrixType& ridgelets, MatrixType& voxels, double lambda);
	void FISTA(MatrixType & x, int N_splits);
	void FISTA2(MatrixType & x, int spl);
};

#endif // !SOLVERS
