#include "SOLVERS.h"

SOLVERS::SOLVERS() {
	cerr << "Minimal set of argumets: ridgelets basis, full DWI array or matrix/vector with voxel(s). "
		"The last parameter - lambda value is optional.\n";
}

SOLVERS::SOLVERS(MatrixType ridgelets, MatrixType voxels) {
	lmd = 0.1;
	A = ridgelets;
	s = voxels;
}

SOLVERS::SOLVERS(MatrixType ridgelets, MatrixType voxels, double lambda) {
	lmd = lambda;
	A = ridgelets;
	s = voxels;
}

MatrixType SOLVERS::FISTA() {
	cout << "Start computing ridgelets coefficients..." << endl;
	MatrixType dummyZeros = MatrixType::Zero(A.cols(), s.cols());
	MatrixType y = dummyZeros;
	MatrixType x_old = dummyZeros;

	double t_old = 1;
	double t = 0;
	double e_old = 1e32;
	double e;

	for (int iter = 0; iter < 10; ++iter) {
		x = y + A.transpose() * (s - A * y);
		//Soft thresholding
		x = ((x.cwiseAbs().array() - lmd).cwiseMax(dummyZeros.array())).cwiseProduct(x.array().sign());

		e = (0.5 * (A * x - s).array().pow(2).colwise().sum().array() +
			lmd * x.cwiseAbs().colwise().sum().array()).maxCoeff();
		//cout << endl << e;
		if ((e_old - e) / e_old < 0.001)
			break;
		else
			e_old = e;

		//Nesterov acceleration
		t = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
		y = x + ((t_old - 1) / t) * (x - x_old);
		x_old = x;
		t_old = t;
		// cout << "Iteration " << iter << " finished\n";
	}
	return x;
}