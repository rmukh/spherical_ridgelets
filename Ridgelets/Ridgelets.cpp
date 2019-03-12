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
	if (res_dmri) {
		return EXIT_SUCCESS;
	}
	MaskImagePointer mask;
	int res_mask = data.readMask("C:\\Users\\mukho\\Desktop\\01009-mask.nhdr", mask);
	if (res_mask) {
		return EXIT_SUCCESS;
	}

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(dMRI, signal, nGradImgs, nOfImgs);

	// just for debugging purpose and development of MatrixType to ITK type image

	DiffusionImagePointer dMRI_new = DiffusionImageType::New();
	data.copy_header(dMRI, dMRI_new);
	dMRI_new->SetNumberOfComponentsPerPixel(signal.rows());
	dMRI_new->Allocate();

	data.Matrix2DWI(dMRI_new, signal);
	data.save_to_file<DiffusionImageType>("C:\\Users\\mukho\\Desktop\\test.nhdr", dMRI_new);

	//cout << "voxel 50000" << signal.col(50000) << endl;

	// Demo data for solvers tests
	//MatrixType g(51, 3);
	//MatrixType s(51, 16);

	//data.readTestData(g, s);

	//icosahedron
	//UtilMath m;

	//MatrixType fcs;
	//MatrixType nu;
	//m.icosahedron(nu, fcs, 1);

	//data.matrixToFile("C:\\Users\\mukho\\Desktop\\fcs.txt", fcs);
	//data.matrixToFile("C:\\Users\\mukho\\Desktop\\u.txt", u);

	//high_resolution_clock::time_point t1 = high_resolution_clock::now(); //start timer point

	//SPH_RIDG ridg(2, 0.5);
	//MatrixType A = ridg.RBasis(g);
	//A = ridg.normBasis(A);

	//SOLVERS slv(A, s, 0.1);
	//MatrixType C = slv.FISTA();
	//cout << endl << C;

	////ODF
	//MatrixType Q = ridg.QBasis(nu); //Build a Q basis
	//MatrixType ODF = C * Q.transpose(); //Computer ODF

	//high_resolution_clock::time_point u2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>(u2 - t1).count();
	//cout << "Execution time " << duration << " seconds" << endl;

	// import values from matlab for tests and debugging
	//MatrixType u;
	//MatrixType fcs;
	//MatrixType W;
	//data.fileToMatrix("C:\\Users\\mukho\\Desktop\\fcs_external.txt", fcs);
	//data.fileToMatrix("C:\\Users\\mukho\\Desktop\\u_external.txt", u);
	//data.fileToMatrix("C:\\Users\\mukho\\Desktop\\odf_one_vol.txt", W);

	//fcs = fcs.array() - 1;

	//vector<vector<unsigned>> conn;
	//m.FindConnectivity(conn, fcs, u.rows());

	//MatrixType ex;
	//MatrixType d;
	//m.FindODFMaxima(ex, d, W, conn, u);
	//cout << d;

	return 0;
}