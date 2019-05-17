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
void SOLVERS<pT, RT, ST>::FISTA(ST& x, int N_splits) {
	cout << "Start computing ridgelets coefficients..." << endl;

	auto start = high_resolution_clock::now();

	x.resize(A->cols(), s->cols());
	unsigned split_size = floor(s->cols() / N_splits);

#pragma omp parallel for
	for (int it = 0; it < N_splits; ++it) {
		ST x_block;
		ST s_block;
		ST y;
		ST x_old;

		s_block = s->block(0, it * split_size, s->rows(), split_size);
		y = ST::Zero(A->cols(), s_block.cols());
		x_old = ST::Zero(A->cols(), s_block.cols());

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

template <class pT, class RT, class ST>
void SOLVERS<pT, RT, ST>::FISTA(ST& x) {
	x.resize(A->cols());
	ST y;
	ST x_old;

	y = ST::Zero(A->cols());
	x_old = ST::Zero(A->cols());

	pT t_old = 1;
	pT t = 0;
	pT e_old = 1e32;
	pT e;

	for (int iter = 0; iter < 2000; ++iter) {
		x = y + A->transpose() * (*s - *A * y);

		//Soft thresholding
		x = ((x.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x.array().sign());

		e = ((0.5 * (*A * x - *s).array().pow(2).colwise().sum().array()) +
			(lmd * x.cwiseAbs().colwise().sum().array())).maxCoeff();

		if ((e_old - e) / e_old < 0.001)
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

template <class pT, class RT, class ST>
void loop_block(ST& x, ST& y, ST& s) {
	pT t_old = 1;
	pT t = 0;
	pT e_old = 1e32;
	pT e;

	for (int iter = 0; iter < 2000; ++iter) {
		x = y + A->transpose() * (*s - *A * y);

		//Soft thresholding
		x = ((x.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x.array().sign());

		e = ((0.5 * (*A * x - *s).array().pow(2).colwise().sum().array()) +
			(lmd * x.cwiseAbs().colwise().sum().array())).maxCoeff();

		if ((e_old - e) / e_old < 0.001)
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