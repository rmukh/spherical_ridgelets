// Ridgelets.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//precompiled header
#include "pch.h"

//user defined classes
#include "UtilMath.h"
#include "SOLVERS.h"
#include "SPH_RIDG.h"
#include "DATA_SOURCE.h"

int main()
{
	Eigen::initParallel();

	DATA_SOURCE data;
	MatrixType GradientDirections(0, 3); // Matrix with dMRI image gradient directions
	DiffusionImagePointer dMRI;
	unsigned nGradImgs = 0; // Number of gradient images
	unsigned nOfImgs = 0; // Total number of images (including b0)
	int res_dmri = data.readVolume("C:\\Users\\mukho\\Desktop\\01009-dwi-Ed.nhdr", GradientDirections, dMRI, nGradImgs, nOfImgs);
	if (res_dmri)
		return EXIT_SUCCESS;

	MaskImagePointer mask;
	int res_mask = data.readMask("C:\\Users\\mukho\\Desktop\\01009-mask.nhdr", mask);
	if (res_mask)
		return EXIT_SUCCESS;

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(dMRI, mask, signal, nGradImgs, nOfImgs);

	UtilMath m;

	MatrixType fcs;
	MatrixType nu;
	m.icosahedron(nu, fcs, 4);

	SPH_RIDG ridg(2, 0.5);
	MatrixType A = ridg.RBasis(GradientDirections);
	A = ridg.normBasis(A);

	SOLVERS slv(A, signal, 0.1);
	MatrixType C = slv.FISTA();

	//ODF
	MatrixType Q = ridg.QBasis(nu); //Build a Q basis
	MatrixType ODF = C * Q.transpose(); //Computer ODF - Might need transpose!!!
	
	// Save to file what user requested through command line
	// R coefficients
	DiffusionImagePointer R_coeff = DiffusionImageType::New();
	data.copy_header(dMRI, R_coeff);
	R_coeff->SetNumberOfComponentsPerPixel(C.rows());
	R_coeff->Allocate();

	data.Matrix2DWI(R_coeff, mask, C);
	data.save_to_file<DiffusionImageType>("C:\\Users\\mukho\\Desktop\\R_coeff.nhdr", R_coeff, false);

	// Q coefficients
	DiffusionImagePointer Q_coeff = DiffusionImageType::New();
	data.copy_header(dMRI, Q_coeff);
	Q_coeff->SetNumberOfComponentsPerPixel(Q.rows());
	Q_coeff->Allocate();

	data.Matrix2DWI(Q_coeff, mask, Q);
	data.save_to_file<DiffusionImageType>("C:\\Users\\mukho\\Desktop\\Q_coeff.nhdr", Q_coeff, false);

	// ODF values
	DiffusionImagePointer ODF_vals = DiffusionImageType::New();
	data.copy_header(dMRI, ODF_vals);
	ODF_vals->SetNumberOfComponentsPerPixel(ODF.rows());
	ODF_vals->Allocate();

	data.Matrix2DWI(ODF_vals, mask, ODF);
	data.save_to_file<DiffusionImageType>("C:\\Users\\mukho\\Desktop\\ODF_vals.nhdr", ODF_vals, false);

	//vector<vector<unsigned>> conn;
	//m.FindConnectivity(conn, fcs, nu.rows());

	//MatrixType ex;
	//MatrixType d;
	//m.FindODFMaxima(ex, d, W, conn, nu);

	return 0;
}