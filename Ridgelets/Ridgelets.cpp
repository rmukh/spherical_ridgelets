// Ridgelets.cpp : This file contains the 'main' function. Program execution begins and ends there.
/*
 Copyright (c) 2019 Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

 * Description:
 *     A C++ implementation of spherical ridgelets and orientation distribution function.
 *     To find spherical ridgelets coefficients A Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) used.
 * Dependencies:
 *     ITK, Eigen
 * Authors:
 *     Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi
*/

//user defined classes
#include "rdgls_types.h"
#include "UtilMath.h"
#include "SOLVERS.h"
#include "SPH_RIDG.h"
#include "DATA_SOURCE.h"

int main(int argc, char* argv[])
{
	Eigen::initParallel();
	DATA_SOURCE data;

	// Parse input parameters from CLI
	DATA_SOURCE::input_parse input_args;
	if (data.CLI(argc, argv, input_args))
		return EXIT_SUCCESS;

	// Set number of threads
	if (input_args.nth != -1) {
		omp_set_num_threads(input_args.nth);
		Eigen::setNbThreads(input_args.nth);
	}

	// Read mask
	MaskImagePointer mask;
	int res_mask = data.readMask(input_args.input_mask, mask);
	if (res_mask)
		return EXIT_SUCCESS;

	// 4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	MatrixType GradientDirections; // Matrix with dMRI image gradient directions
	int res_dmri = data.DWI2Matrix(input_args.input_dmri, mask, signal, GradientDirections);
	if (res_dmri)
		return EXIT_SUCCESS;

	// Beginning of the main computational part
	SPH_RIDG<precisionType, MatrixType, VectorType>ridg(input_args.sph_J, 1/input_args.sph_rho);
	MatrixType A;
	ridg.RBasis(A, GradientDirections);
	ridg.normBasis(A);

	if (input_args.n_splits == -1)
		input_args.n_splits = data.compute_splits(signal.cols());
	
	data.estimate_memory(signal, A, input_args);
	data.short_summary(input_args);

	MatrixType C;
	{
		SOLVERS<precisionType, MatrixType, MatrixType> slv(A, signal, input_args.fista_lambda);
		slv.FISTA(C, input_args.n_splits);  //have a potentinal for optimization
	}

	// Save to file(s) user requested through command line
	// Ridgelets coefficients
	if (!input_args.output_ridgelets.empty()) {
		cout << "Saving ridgelets coefficients..." << endl;
		DiffusionImagePointer Ridg_coeff = DiffusionImageType::New();
		data.set_header(Ridg_coeff);
		Ridg_coeff->SetNumberOfComponentsPerPixel(C.rows());
		Ridg_coeff->Allocate();

		data.Matrix2DWI(Ridg_coeff, mask, C);
		data.save_to_file<DiffusionImageType>(input_args.output_ridgelets, Ridg_coeff, input_args.is_compress);
	}

	// A*c (signal recon)
	if (!input_args.signal_recon.empty()) {
		cout << "Saving signal reconstruction..." << endl;
		MatrixType SR = A * C;

		DiffusionImagePointer s_coeff = DiffusionImageType::New();
		data.set_header(s_coeff);
		s_coeff->SetNumberOfComponentsPerPixel(SR.rows());
		s_coeff->Allocate();

		data.Matrix2DWI(s_coeff, mask, SR);
		data.save_to_file<DiffusionImageType>(input_args.signal_recon, s_coeff, input_args.is_compress);
	}

	UtilMath<precisionType, MatrixType, VectorType> m;
	MatrixType fcs;
	MatrixType nu;
	MatrixType Q;

	if (!input_args.output_odf.empty() || !input_args.output_fiber_max_odf.empty()) {
		m.icosahedron(nu, fcs, input_args.lvl);
		ridg.QBasis(Q, nu); //Build a Q basis
	}

	// ODF volume
	if (!input_args.output_odf.empty()) {
		MatrixType ODF = Q * C;
		cout << "Saving ODF values..." << endl;
		DiffusionImagePointer ODF_vals = DiffusionImageType::New();
		data.set_header(ODF_vals);
		ODF_vals->SetNumberOfComponentsPerPixel(Q.rows());
		ODF_vals->Allocate();

		data.Matrix2DWI(ODF_vals, mask, ODF);
		data.save_to_file<DiffusionImageType>(input_args.output_odf, ODF_vals, input_args.is_compress);
	}

	// Maximum directions and values of ODF
	if (!input_args.output_fiber_max_odf.empty()) {
		MatrixType ex_d;
		vector<vector<unsigned>> conn;

		m.FindConnectivity(conn, fcs, nu.rows());
		m.FindMaxODFMaxInDMRI(ex_d, Q, C, conn, nu, input_args.max_odf_thresh);

		cout << "Saving maxima ODF direction and value..." << endl;
		DiffusionImagePointer modf = DiffusionImageType::New();
		data.set_header(modf);
		modf->SetNumberOfComponentsPerPixel(ex_d.rows());
		modf->Allocate();

		data.Matrix2DWI(modf, mask, ex_d);
		data.save_to_file<DiffusionImageType>(input_args.output_fiber_max_odf, modf, input_args.is_compress);
	}

	return EXIT_SUCCESS;
}
