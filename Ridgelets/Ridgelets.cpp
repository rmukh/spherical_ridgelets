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
			cout << U << endl;
			/*
				Find unique rows 
				Not the fanciest and most optimal way, but
				1) icosahedron function does not intend to be called many times
				2) uses kind of hashtable as serious boys usually do :)
			*/

			// Define hashtable
			unordered_map<string, bool> hTable;

			// Define array of existing indicies
			vector<int> duplicates;

			// Preallocate string for faster string concatenation
			string key;
			size_t added_length = 3 * to_string(U(0, 0)).length();
			key.reserve(key.length() + added_length);

			// Iterate over matrix
			for (unsigned i = 0; i < U.rows(); ++i) {
				// Create unique key from row valuse
				key = to_string(U(i, 0)).append(to_string(U(i, 1))).append(to_string(U(i, 2))); 
				
				// If element exists in hash table
				if (hTable.count(key) > 0) 
					duplicates.push_back(i);
				else
					hTable.insert(pair<string, bool>(key, true));
			}

			for (const int i : duplicates)
				cout << i << ' ';
			cout << endl;

			Map<VectorXi> rowsToKeep(duplicates.data(), duplicates.size());
			VectorXd allCols = VectorXd::LinSpaced(3, 0, 2);
			cout << U.outerStride();
			Map<MatrixType, 0, OuterStride<>> U2(U.data(), rowsToKeep.rows(), 3, OuterStride<>(U.outerStride() * 3));
		}
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
