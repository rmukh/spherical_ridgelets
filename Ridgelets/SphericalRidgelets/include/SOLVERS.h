#ifndef SOLVERS_H
#define SOLVERS_H

#include "rdgls_types.h"

//precisionType - pT, MatrixType - MT, VectorType - VT
template <class pT, class RT, class ST>
class SOLVERS
{
public:
	RT *A; //rigelets basis
	ST *s; //input voxel(s)
	pT lmd; //labmda parameter for FISTA training

	SOLVERS();
	~SOLVERS();

	SOLVERS(RT& ridgelets, ST& voxels);
	SOLVERS(RT& ridgelets, ST& voxels, pT lambda);

	void FISTA(ST& x, int N_splits, int n_iterations, precisionType tolerance);
	void FISTA(ST& x, int n_iterations, precisionType tolerance);

	void loop_block(ST& x, ST& sig, int n_iterations, precisionType tolerance);
	ST SolveSingle_py(RT& basis, ST& sig, precisionType lmd, int n_iterations, precisionType tolerance);
};

#include "SOLVERS.hpp"

#endif // !SOLVERS
