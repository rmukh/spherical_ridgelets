#include "SOLVERS.h"

#ifndef SOLVERS_IMPL
#define SOLVERS_IMPL

template <class pT, class MT>
SOLVERS<pT, MT>::SOLVERS() : A(NULL), s(NULL) {
	cerr << "Minimal set of argumets: ridgelets basis, full DWI array or matrix/vector with voxel(s). "
		"The last parameter - lambda value is optional.\n";
}

template <class pT, class MT>
SOLVERS<pT, MT>::~SOLVERS() {}

template <class pT, class MT>
SOLVERS<pT, MT>::SOLVERS(MT& ridgelets, MT& voxels) : A(&ridgelets), s(&voxels), lmd(0.1) {}

template <class pT, class MT>
SOLVERS<pT, MT>::SOLVERS(MT& ridgelets, MT& voxels, pT lambda) : A(&ridgelets), s(&voxels), lmd(lambda) {}

template <class pT, class MT>
void SOLVERS<pT, MT>::FISTA(MT& x, int N_splits) {
	cout << "Start computing ridgelets coefficients..." << endl;

	auto start = high_resolution_clock::now();

	x.resize(A->cols(), s->cols());
	unsigned split_size = floor(s->cols() / N_splits);

	#pragma omp parallel for
	for (int it = 0; it < N_splits; ++it) {
		MT x_block;
		MT s_block;
		MT y;
		MT x_old;

		s_block = s->block(0, it * split_size, s->rows(), split_size);
		y = MT::Zero(A->cols(), s_block.cols());
		x_old = MT::Zero(A->cols(), s_block.cols());

		pT t_old = 1;
		pT t = 0;
		pT e_old = 1e32;
		pT e;

		for (int iter = 0; iter < 2000; ++iter) {
			x_block = y + A->transpose() * (s_block - *A * y);

			//Soft thresholding
			x_block = ((x_block.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x_block.array().sign());

			e = ((0.5 * (*A * x_block - s_block).array().pow(2).colwise().sum().array()) +
				(lmd * x_block.cwiseAbs().colwise().sum().array())).maxCoeff();

			if ((e_old - e) / e_old < 0.001)
				break;
			else
				e_old = e;

			//Nesterov acceleration
			t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
			y = x_block + ((t_old - 1) / t) * (x_block - x_old);
			x_old = x_block;
			t_old = t;
		}
		x.block(0, it * split_size, x_block.rows(), x_block.cols()) = x_block;
	}

	auto stop = high_resolution_clock::now();
	auto ds = duration_cast<seconds>(stop - start);
	auto dm = duration_cast<minutes>(stop - start);
	cout << "Computations took " << ds.count() << " seconds ~ " << dm.count() << "+ minutes" << endl;
}

#endif