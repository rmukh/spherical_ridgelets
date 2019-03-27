#include "SOLVERS.h"

SOLVERS::SOLVERS() : A(nullptr), s(nullptr) {
	cerr << "Minimal set of argumets: ridgelets basis, full DWI array or matrix/vector with voxel(s). "
		"The last parameter - lambda value is optional.\n";
}

SOLVERS::~SOLVERS() { cout << "SOLVER destructed" << endl; }

SOLVERS::SOLVERS(MatrixType& ridgelets, MatrixType& voxels) : A(&ridgelets), s(&voxels), lmd(0.1) {}

SOLVERS::SOLVERS(MatrixType& ridgelets, MatrixType& voxels, double lambda) : A(&ridgelets), s(&voxels), lmd(lambda) {}

void SOLVERS::FISTA(MatrixType& x, int N_splits) {
	cout << "Start computing ridgelets coefficients..." << endl;

	N_splits = s->cols() / 2;
	x.resize(A->cols(), s->cols());
	unsigned split_size = floor(s->cols() / N_splits);

	#pragma omp parallel for
	for (int it = 0; it < N_splits; ++it) {
		MatrixType x_block;
		MatrixType s_block;
		MatrixType y;
		MatrixType x_old;

		s_block = s->block(0, it * split_size, s->rows(), split_size);
		y = MatrixType::Zero(A->cols(), s_block.cols());
		x_old = MatrixType::Zero(A->cols(), s_block.cols());

		double t_old = 1;
		double t = 0;
		double e_old = 1e32;
		double e;

		for (int iter = 0; iter < 2000; ++iter) {
			x_block = y + A->transpose() * (s_block - *A * y);

			//Soft thresholding
			x_block = ((x_block.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x_block.array().sign());

			e = ((0.5 * (*A * x_block - s_block).array().pow(2).colwise().sum().array()) +
				(lmd * x_block.cwiseAbs().colwise().sum().array())).maxCoeff();
			//cout << e << endl;
			if ((e_old - e) / e_old < 0.001)
				break;
			else
				e_old = e;

			//Nesterov acceleration
			t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
			y = x_block + ((t_old - 1) / t) * (x_block - x_old);
			x_old = x_block;
			t_old = t;
			//cout << "Iteration " << iter << " in split " << it << endl;
		}
		x.block(0, it * split_size, x_block.rows(), x_block.cols()) = x_block;
	}
}

void SOLVERS::FISTA2(MatrixType& x, int spl) {
	cout << "Start computing ridgelets coefficients..." << endl;

	MatrixType y = MatrixType::Zero(A->cols(), s->cols());
	MatrixType x_old = MatrixType::Zero(A->cols(), s->cols());

	double t_old = 1;
	double t = 0;
	double e_old = 1e32;
	double e;

	for (int iter = 0; iter < 2000; ++iter) {
		x = y + A->transpose() * (*s - *A * y);

		//Soft thresholding
		x = ((x.cwiseAbs().array() - lmd).cwiseMax(0)).cwiseProduct(x.array().sign());

		e = ((0.5 * (*A * x - *s).array().pow(2).colwise().sum().array()) +
			(lmd * x.cwiseAbs().colwise().sum().array())).maxCoeff();
		//cout << "e " << e_old << " " << e << " " << (e_old - e) / e_old << endl;
		if ((e_old - e) / e_old < 0.001)
			break;
		else
			e_old = e;

		//Nesterov acceleration
		t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
		y = x + ((t_old - 1) / t) * (x - x_old);
		x_old = x;
		t_old = t;
		//cout << "Iteration " << iter << " finished" << endl;
	}
}
