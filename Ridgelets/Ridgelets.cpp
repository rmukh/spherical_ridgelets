// Ridgelets.cpp : This file contains the 'main' function. Program execution begins and ends there.
/*
 Copyright (c) 2021 Rinat Mukhometzianov, Oleg Michailovich, Yogesh Rathi

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
#if defined(USE_FLOAT) || defined(UKF_USE_FLOAT)
	std::cout << "The package was compiled for float precision computations" << std::endl;
#endif

	Eigen::initParallel();

#if defined(_OPENMP)
	cout << "OpenMP enabled" << endl;
#else
	cout << "No OpenMP!" << endl;
#endif

	DATA_SOURCE data;

	// Parse input parameters from CLI
	DATA_SOURCE::input_parse input_args {
		"","","","","","","","","",
		0.7,0.01,3.125,4,2,-1,-1,false,2000,0.001
	};

	if (data.CLI(argc, argv, &input_args))
		return EXIT_SUCCESS;

#if defined(_OPENMP)
	// Set number of threads
	if (input_args.nth != -1)
	{
		omp_set_num_threads(input_args.nth);
		Eigen::setNbThreads(input_args.nth);
	}
#endif
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
	SPH_RIDG<precisionType, MatrixType, VectorType> ridg(input_args.sph_J, 1 / input_args.sph_rho);
	MatrixType A;
	cout << "Start computing R basis..." << endl;
	ridg.RBasis(A, GradientDirections);
	ridg.normBasis(A);
	
	// Read external gradients
	MatrixType ext_g;
	MatrixType ext_A;
	if (!input_args.external_gradients.empty()) {
		data.fileGradientsToMatrix(input_args.external_gradients, ext_g);
		ridg.RBasis(ext_A, ext_g);
		ridg.normBasis(ext_A);
	}

	if (!input_args.output_A.empty()) // -A
	{
		cout << "Saving A basis..." << endl;
		data.matrixToFile(input_args.output_A, A);
	}

	// If only A output required
	if (!input_args.output_A.empty() &&
		input_args.output_ridgelets.empty() &&
		input_args.signal_recon.empty() &&
		input_args.output_odf.empty() &&
		input_args.output_fiber_max_odf.empty())
	{
		return EXIT_SUCCESS;
	}

	if (input_args.n_splits == -1)
		input_args.n_splits = data.compute_splits(signal.cols());

	data.estimate_memory(signal, A, input_args);
	data.short_summary(input_args);

	MatrixType C;
	{
		SOLVERS<precisionType, MatrixType, MatrixType> slv(A, signal, input_args.fista_lambda);
		cout << "Start computing ridgelets coefficients..." << endl;
		auto start = high_resolution_clock::now();
		//have a potentinal for optimization
		slv.FISTA(C, input_args.n_splits, input_args.fista_iterations, input_args.fista_tolerance); 
		auto stop = high_resolution_clock::now();
		auto ds = duration_cast<seconds>(stop - start);
		auto dm = duration_cast<minutes>(stop - start);
		cout << "Computations took " << ds.count() << " seconds ~ " << dm.count() << "+ minutes" << endl;
	}

	// Save to file(s) user requested through command line
	// Ridgelets coefficients
	if (!input_args.output_ridgelets.empty()) // -ridg
	{
		cout << "Saving ridgelets coefficients..." << endl;
		DiffusionImagePointer Ridg_coeff = DiffusionImageType::New();
		data.set_header(Ridg_coeff);
		Ridg_coeff->SetNumberOfComponentsPerPixel(C.rows());
		Ridg_coeff->Allocate();

		data.Matrix2DWI(Ridg_coeff, mask, C);
		data.save_to_file<DiffusionImageType>(input_args.output_ridgelets, Ridg_coeff, input_args.is_compress);
	}

	// A*c (signal recon)
	if (!input_args.signal_recon.empty()) // -sr
	{
		cout << "Saving signal reconstruction..." << endl;
		MatrixType SR = A * C;

		DiffusionImagePointer s_coeff = DiffusionImageType::New();
		data.set_header(s_coeff);
		s_coeff->SetNumberOfComponentsPerPixel(SR.rows());
		s_coeff->Allocate();

		data.Matrix2DWI(s_coeff, mask, SR);
		data.save_to_file<DiffusionImageType>(input_args.signal_recon, s_coeff, input_args.is_compress);
	}

	// ext_A*c (signal recon with an external gradients)
	if (!input_args.ext_signal_recon.empty()) // -ext_sr
	{
		cout << "Saving signal reconstruction with an external gradients..." << endl;
		MatrixType SR = ext_A * C;

		DiffusionImagePointer s_coeff = DiffusionImageType::New();
		data.set_header(s_coeff);
		s_coeff->SetNumberOfComponentsPerPixel(SR.rows());
		s_coeff->Allocate();

		data.Matrix2DWI(s_coeff, mask, SR);
		data.save_to_file<DiffusionImageType>(input_args.ext_signal_recon, s_coeff, input_args.is_compress);
	}

	UtilMath<precisionType, MatrixType, VectorType> m;
	MatrixType fcs;
	MatrixType nu;
	MatrixType Q;

	if (!input_args.output_odf.empty() || !input_args.output_fiber_max_odf.empty())
	{
		cout << "Start computing icosahedron..." << endl;
		m.icosahedron(nu, fcs, input_args.lvl);
		cout << "Start computing Q basis..." << endl;
		ridg.QBasis(Q, nu); //Build a Q basis
	}

	// ODF volume
	if (!input_args.output_odf.empty()) // -odf
	{
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
	if (!input_args.output_fiber_max_odf.empty()) // -omd
	{
		cout << "Start computing ODF maxima..." << endl;
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
