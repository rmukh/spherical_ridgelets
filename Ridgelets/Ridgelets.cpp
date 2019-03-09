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
	DATA_SOURCE data;
	/*
	MatrixType GradientDirections(0, 3);
	DiffusionImagePointer image;
	unsigned nGradImgs = 0;
	unsigned nOfImgs = 0;
	data.readNRRD("C:\\Users\\renat\\Desktop\\01009-dwi-Ed.nhdr", GradientDirections, image, nGradImgs, nOfImgs);

	const DiffusionImageType::IndexType pixelIndex = { {27,29,37} }; //Position {X,Y,Z}
	DiffusionImageType::PixelType value = image->GetPixel(pixelIndex);

	cout << "Gradient Directions: \n" << GradientDirections << endl;
	cout << value << endl;

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(image, signal, nGradImgs, nOfImgs);
	image = nullptr;

	//cout << "voxel 50000" << signal.col(50000) << endl;
	*/
	MatrixType g(51, 3);
	MatrixType s(51, 16);

	data.readTestData(g, s);

	// import values from matlab
	MatrixType u;
	MatrixType fcs;
	data.fileToMatrix("C:\\Users\\mukho\\Desktop\\fcs_external.txt", fcs);
	data.fileToMatrix("C:\\Users\\mukho\\Desktop\\u_external.txt", u);

	//icosahedron
	UtilMath m;

	//MatrixType fcs;
	//MatrixType u;
	//m.icosahedron(u, fcs, 1);

	//data.matrixToFile("C:\\Users\\mukho\\Desktop\\fcs.txt", fcs);
	//data.matrixToFile("C:\\Users\\mukho\\Desktop\\u.txt", u);

	//for (vtkIdType i = 0; i < raw->GetNumberOfPoints(); i++)
	//{
	//	double p[3];
	//	raw->GetCell(i, p);
	//	cout << "point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << endl;
	//}

	//
	//high_resolution_clock::time_point t1 = high_resolution_clock::now(); //start timer point

	//SPH_RIDG ridg(2, 0.5);
	//MatrixType A = ridg.RBasis(g);
	//A = ridg.normBasis(A)

	//SOLVERS slv(A, s, 0.1);
	//MatrixType C = slv.FISTA();
	//cout << endl << C;

	////ODF
	////MatrixType Q = ridg.QBasis(nu); //Build a Q basis
	////MatrixType ODF = C * Q.transpose(); //Computer ODF

	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>(t2 - t1).count();
	//cout << "Execution time " << duration << " seconds" << endl;

	// FindConnectivity
	unsigned N = u.rows();

	vector<Eigen::Index> a1;
	m.column_find(a1, fcs, 0, true, 0+1);

	vector<Eigen::Index> a2;
	m.column_find(a2, fcs, 1, true, 0+1);

	vector<Eigen::Index> a3;
	m.column_find(a3, fcs, 2, true, 0+1);

	data.printVec("a1", a1);
	data.printVec("a2", a2);
	data.printVec("a3", a3);

	MatrixType t(a3.size(), 2);
	for (unsigned i = 0; i < a3.size(); ++i)
		t.row(i) = fcs.block<1, 2>(a3.at(i), 0);

	t.resize(2 * a3.size(), 1);
	vector<int> un;
	m.unique_sorted(un, t);

	cout << "u " << endl << t << endl;
	data.printVec("uniques", un);

	return 0;
}
