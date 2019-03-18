// Ridgelets.cpp : This file contains the 'main' function. Program execution begins and ends there.

//user defined classes
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

	MatrixType GradientDirections(0, 3); // Matrix with dMRI image gradient directions
	DiffusionImagePointer dMRI;
	unsigned nGradImgs = 0; // Number of gradient images
	unsigned nOfImgs = 0; // Total number of images (including b0)
	int res_dmri = data.readVolume(input_args.input_dmri, GradientDirections, dMRI, nGradImgs, nOfImgs);
	if (res_dmri)
		return EXIT_SUCCESS;

	MaskImagePointer mask;
	int res_mask = data.readMask(input_args.input_mask, mask);
	if (res_mask)
		return EXIT_SUCCESS;

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(dMRI, mask, signal, nGradImgs, nOfImgs);

	// Beginning of the main computational part
	SPH_RIDG ridg(2, 0.5);
	MatrixType A = ridg.RBasis(GradientDirections);
	A = ridg.normBasis(A);

	SOLVERS slv(A, signal, 0.1);
	MatrixType C = slv.FISTA();

	// Save to file what user requested through command line
	// Ridgelets coefficients
	if (!input_args.output_ridgelets.empty()) {
		cout << "Saving ridgelets coefficients..." << endl;
		DiffusionImagePointer Ridg_coeff = DiffusionImageType::New();
		data.copy_header(dMRI, Ridg_coeff);
		Ridg_coeff->SetNumberOfComponentsPerPixel(C.rows());
		Ridg_coeff->Allocate();

		data.Matrix2DWI(Ridg_coeff, mask, C);
		data.save_to_file<DiffusionImageType>(input_args.output_ridgelets, Ridg_coeff, input_args.is_compress);
	}

	UtilMath m;
	MatrixType fcs;
	MatrixType nu;
	MatrixType Q;
	MatrixType ODF;

	if (!input_args.output_odf.empty() || !input_args.output_fiber_max_odf.empty() || !input_args.output_dirs.empty()) {
		m.icosahedron(nu, fcs, 1);
		Q = ridg.QBasis(nu); //Build a Q basis
		ODF = Q * C;
	}

	// ODF volume
	if (!input_args.output_odf.empty()) {
		cout << "Saving ODF values..." << endl;
		DiffusionImagePointer ODF_vals = DiffusionImageType::New();
		data.copy_header(dMRI, ODF_vals);
		ODF_vals->SetNumberOfComponentsPerPixel(Q.rows());
		ODF_vals->Allocate();

		data.Matrix2DWI(ODF_vals, mask, ODF);
		data.save_to_file<DiffusionImageType>(input_args.output_odf, ODF_vals, input_args.is_compress);
	}

	MatrixType ex;
	MatrixType d;
	vector<vector<unsigned>> conn;

	if (!input_args.output_fiber_max_odf.empty() || !input_args.output_dirs.empty()) {
		m.FindConnectivity(conn, fcs, nu.rows());
		m.FindMaxODFMaxInDMRI(ex, d, ODF, conn, nu);
	}

	if (!input_args.output_fiber_max_odf.empty()) {
		cout << "Saving maxima ODF..." << endl;
		MatrixType max_odf = MatrixType::Zero(ex.rows(), ex.cols());
		for (unsigned i = 0; i < ex.cols(); ++i)
			for (unsigned j = 0; j < ex.rows(); ++j)
				max_odf(j, i) = ODF(ex(j, i), i);

		DiffusionImagePointer modf = DiffusionImageType::New();
		data.copy_header(dMRI, modf);
		modf->SetNumberOfComponentsPerPixel(max_odf.rows());
		modf->Allocate();

		data.Matrix2DWI(modf, mask, max_odf);
		data.save_to_file<DiffusionImageType>(input_args.output_odf, modf, input_args.is_compress);
	}

	if (!input_args.output_dirs.empty()) {
		cout << "Saving directions..." << endl;
		DiffusionImagePointer dirs = DiffusionImageType::New();
		data.copy_header(dMRI, dirs);
		dirs->SetNumberOfComponentsPerPixel(d.rows());
		dirs->Allocate();

		data.Matrix2DWI(dirs, mask, d);
		data.save_to_file<DiffusionImageType>(input_args.output_odf, dirs, input_args.is_compress);
	}
	return EXIT_SUCCESS;
}