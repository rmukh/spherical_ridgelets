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

	//icosahedron
	UtilMath m;
	unsigned level = 1;

	double C = 1 / sqrt(1.25);
	MatrixType t = (2 * m.PI / 5.0) * VectorXd::LinSpaced(5, 0, 4);
	MatrixType u1(5, 3);
	u1 << C * t.array().cos(), C * t.array().sin(), C * 0.5 * MatrixType::Ones(5, 1);
	MatrixType u2(5, 3);
	u2 << C * (t.array() + 0.2 * m.PI).cos(), C * (t.array() + 0.2 * m.PI).sin(), -0.5 * C * MatrixType::Ones(5, 1);
	MatrixType u(12, 3);
	u << 0, 0, 1, u1, u2, 0, 0, -1;
	cout << endl;
	/* */
	if (level > 0) {
		for (unsigned lev = 1; lev <= level; ++lev) {
			MatrixType fcs = m.convhulln(u);
			unsigned N = fcs.rows();
			MatrixType U = MatrixType::Zero(3 * N, 3);
			MatrixType A;
			MatrixType B;
			MatrixType C;
			for (unsigned k = 0; k < N; ++k) {
				A = u.row(fcs(k, 0));
				B = u.row(fcs(k, 1));
				C = u.row(fcs(k, 2));
				U.block<3, 3>(3 * k, 0) << 0.5 * (A + B), 0.5 * (B + C), 0.5 * (A + C);
			}

			vector<int> uniques;
			m.unique_rows(uniques, U);

			// Normalize and add to u
			unsigned u_len = u.rows();
			unsigned unique_len = uniques.size();
			u.conservativeResize(u.rows() + unique_len, u.cols());

			for (unsigned i = 0, j = 0 + u_len; i < unique_len, j < unique_len + u_len; ++i, ++j)
				u.row(j) = U.row(uniques.at(i)) / U.row(uniques.at(i)).norm();
		}
		cout << u << endl << endl;

		// Matrix column to std vector
		vector<double> uc3;
		uc3.resize(u.col(2).size());
		VectorXd::Map(&uc3[0], u.col(2).size()) = u.col(2);

		// Mapping from value to index
		multimap<double, unsigned> indx;
		for (auto it = uc3.rbegin(); it != uc3.rend(); ++it)
			indx.insert(make_pair(*it, it - uc3.rbegin())); //find shift by 1 element

		for (auto it = indx.begin(); it != indx.end(); ++it)
			cout << it->second << endl;

		for (unsigned i = 0; i < uc3.size(); ++i)
			cout << uc3.at(i) << " ";
	}

	//MatrixType fcs = m.convhulln(u);
	//cout << "faces" << endl << fcs.array() + 1;

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

	return 0;
}
