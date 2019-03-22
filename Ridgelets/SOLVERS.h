#ifndef SOLVERS_H
#define SOLVERS_H

#include "UtilMath.h"

class SOLVERS
{
public:
	MatrixType x; //output vector/matrix/4D array of ridgelets coefficients
	MatrixType *A; //rigelets basis
	MatrixType *s; //input voxel(s)
	double lmd; //labmda parameter for FISTA training

	SOLVERS();
	~SOLVERS();

	SOLVERS(MatrixType& ridgelets, MatrixType& voxels);
	SOLVERS(MatrixType& ridgelets, MatrixType& voxels, double lambda);
	MatrixType FISTA();
};

#endif // !SOLVERS
