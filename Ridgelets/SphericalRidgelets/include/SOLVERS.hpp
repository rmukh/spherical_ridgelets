#include "SOLVERS.h"

#ifndef SOLVERS_IMPL
#define SOLVERS_IMPL

template <class pT, class RT, class ST>
SOLVERS<pT, RT, ST>::SOLVERS() : A(NULL), s(NULL) {
	cerr << "Minimal set of argumets: ridgelets basis, full DWI array or matrix/vector with voxel(s). "
		"The last parameter - lambda value is optional.\n";
}

template <class pT, class RT, class ST>
SOLVERS<pT, RT, ST>::~SOLVERS() {}

// For Matrix (multiple voxels)
template <class pT, class RT, class ST>
SOLVERS<pT, RT, ST>::SOLVERS(RT& ridgelets, ST& voxels) : A(&ridgelets), s(&voxels), lmd(0.1) {}

template <class pT, class RT, class ST>
SOLVERS<pT, RT, ST>::SOLVERS(RT& ridgelets, ST& voxels, pT lambda) : A(&ridgelets), s(&voxels), lmd(lambda) {}

template <class pT, class RT, class ST>
void SOLVERS<pT, RT, ST>::FISTA(ST& x, int N_splits, int n_iterations, precisionType tolerance) {
	x.resize(A->cols(), s->cols());
	int split_size = floor(s->cols() / N_splits);

#if defined(_OPENMP)
	#pragma omp parallel for
#endif
	for (int it = 0; it < N_splits; ++it) {
		ST x_block;
		ST s_block;

		s_block = s->block(0, it * split_size, s->rows(), split_size);
		
		loop_block(x_block, s_block, n_iterations, tolerance);
		x.block(0, it * split_size, x_block.rows(), x_block.cols()) = x_block;
	}
}

template <class pT, class RT, class ST>
void SOLVERS<pT, RT, ST>::FISTA(ST& x, int n_iterations, precisionType tolerance) {
	if (x.cols() == 1)
		x.resize(A->cols(), 1);
	else
		x.resize(A->cols(), s->cols());
	loop_block(x, *s, n_iterations, tolerance);
}

template<class pT, class RT, class ST>
void SOLVERS<pT, RT, ST>::loop_block(ST & x, ST & sig, int n_iterations, precisionType tolerance)
{
	ST y;
	y = ST::Zero(A->cols(), sig.cols());

	ST x_old;
	x_old = ST::Zero(A->cols(), sig.cols());

	pT t_old = 1;
	pT t = 0;
	pT e_old = 1e32;
	pT e;

	for (int iter = 0; iter < n_iterations; ++iter) {
		x = y + A->transpose() * (sig - *A * y);

		//Soft thresholding
		x = ((x.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x.array().sign());

		e = ((0.5 * (*A * x - sig).array().pow(2).colwise().sum().array()) +
			(lmd * x.cwiseAbs().colwise().sum().array())).maxCoeff();

		if (std::fabs(e_old - e) / e_old < tolerance)
			break;
		else
			e_old = e;
		
		//Nesterov acceleration
		t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
		y = x + ((t_old - 1) / t) * (x - x_old);
		x_old = x;
		t_old = t;
	}
}

#endif