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

	// Estimate dMRI Eigen matrix size
	DiffusionImageType::SizeType dmri_size = dMRI->GetLargestPossibleRegion().GetSize();
	long unsigned int dmri_memory = dmri_size[0] * dmri_size[1] * dmri_size[2] 
		* dMRI->GetNumberOfComponentsPerPixel() * sizeof(double);
	
	// Estimate memory consumption by FISTA solver

	if (!input_args.output_odf.empty() || !input_args.output_fiber_max_odf.empty()) {
	}

	cout << "To successfully finish computations you need at least " << dmri_memory / pow(1024, 3) << " GB of RAM and virtual memory combined" << endl;

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(dMRI, mask, signal, nGradImgs, nOfImgs);
	dMRI = nullptr;

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

	if (!input_args.output_odf.empty() || !input_args.output_fiber_max_odf.empty()) {
		m.icosahedron(nu, fcs, input_args.lvl);
		Q = ridg.QBasis(nu); //Build a Q basis
		ODF = Q * C;
	}
	cout << "odf output path " << input_args.output_fiber_max_odf << endl;

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
	
	if (!input_args.output_fiber_max_odf.empty()) {
		MatrixType ex_d;
		vector<vector<unsigned>> conn;

		MatrixType c;

		m.FindConnectivity(conn, fcs, nu.rows());
		m.FindMaxODFMaxInDMRI(ex_d, c, ODF, conn, nu);

		cout << "Saving maxima ODF direction and value..." << endl;
		DiffusionImagePointer modf = DiffusionImageType::New();
		data.copy_header(dMRI, modf);
		modf->SetNumberOfComponentsPerPixel(ex_d.rows());
		modf->Allocate();

		data.Matrix2DWI(modf, mask, ex_d);
		data.save_to_file<DiffusionImageType>(input_args.output_fiber_max_odf, modf, input_args.is_compress);

		//
		DiffusionImagePointer co = DiffusionImageType::New();
		data.copy_header(dMRI, co);
		co->SetNumberOfComponentsPerPixel(c.rows());
		co->Allocate();

		data.Matrix2DWI(co, mask, c);
		data.save_to_file<DiffusionImageType>("C:\\Users\\mukho\\Desktop\\count.nrrd", co, input_args.is_compress);
	}

	return EXIT_SUCCESS;
}
