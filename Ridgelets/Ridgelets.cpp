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
	MatrixType W;
	data.fileToMatrix("C:\\Users\\mukho\\Desktop\\fcs_external.txt", fcs);
	data.fileToMatrix("C:\\Users\\mukho\\Desktop\\u_external.txt", u);
	data.fileToMatrix("C:\\Users\\mukho\\Desktop\\odf_one_vol.txt", W);

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

	//high_resolution_clock::time_point u2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>(u2 - t1).count();
	//cout << "Execution time " << duration << " seconds" << endl;

	vector<vector<unsigned>> conn;
	m.FindConnectivity(conn, fcs, u.rows());

	//FindODFMaxima
	float thresh = 0.7;

	// Standart min-max normalization
	double W_min = W.minCoeff();
	W = (W.array() - W_min) / (W.maxCoeff() - W_min);

	// Find maxima above this point
	MatrixType used = MatrixType::Zero(W.rows(), W.cols());

	// used(W <= thresh) = 1
	for (unsigned i = 0; i < W.size(); ++i)
		if (W(i) <= thresh)
			used(i) = 1;

	unsigned ct = 0;

	MatrixType extrema(0, 0);
	//extrema(0) = 1;

	for (unsigned n = 0; n < used.size(); ++n) {
		if (used(n) == 0) {
			unsigned j = n;
			bool reached_maxima = false;
			while (!reached_maxima) {
				// if (any(W(conn(j).elem) >= W(j)))
				bool if_any = false;
				unsigned conn_row_length = conn[j].size();
				for (unsigned i = 0; i < conn_row_length; i++) {
					if (W(conn[j][i] - 1) >= W(j)) { //remove -1 for real data
						if_any = true;
						break; // trick to speed up computations.
					}
				}
				if (if_any) {
					// [maxw id] = max(W(conn(j).elem))
					unsigned id = 0;
					unsigned maxw = 0;
					for (unsigned i = 0; i < conn_row_length; i++) {
						if (W(conn[j][i] - 1) > maxw) { //remove -1 for real data
							maxw = W(conn[j][i] - 1);
							id = i;
						}
					}

					// We have already traveled this path
					if (used(conn[j][id] - 1)) //remove -1 for real data
						reached_maxima = true;

					// used(conn(j).elem) = 1
					for (unsigned i = 0; i < conn_row_length; ++i)
						used(conn[j][i] - 1) = 1; //remove -1 for real data

					// j = conn(j).elem(id)
					j = conn[j][id] - 1; //remove -1 for real data
				}
				else {
					reached_maxima = true;
					extrema.conservativeResize(1, extrema.cols() + 1);
					extrema(ct) = j;
					for (unsigned i = 0; i < conn_row_length; ++i)
						used(conn[j][i] - 1) = 1; //remove -1 for real data
					ct += 1;
				}
			}
		}
	}

	vector<unsigned> u_extrema;
	m.unique_sorted(u_extrema, extrema);

	unsigned u_extrema_length = u_extrema.size();
	MatrixType directions(u_extrema_length, 3);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		directions.row(i) = u.row(u_extrema.at(i));

	//sort(W(extrema),'descend')
	MatrixType W_e(u_extrema_length, 1);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		W_e(i) = W(u_extrema.at(i));

	multimap<double, unsigned> idxies;
	m.ind_sort_vec(W_e, idxies);

	// Go in reverse order to get descend order as needed in that part of code
	unsigned i = 0;
	vector<unsigned> idx;
	MatrixType directions_sorted(directions.rows(), directions.cols());
	for (multimap<double, unsigned>::reverse_iterator it = idxies.rbegin(); it != idxies.rend(); ++it) {
		directions_sorted.row(i) = directions.row(it->second);
		idx.push_back(it->second);
		i += 1;
	}

	//reverse idx to make in right order
	reverse(idx.begin(), idx.end());

	MatrixType d = MatrixType::Zero(unsigned(extrema.size() / 2) * 2, 3);
	MatrixType ex = MatrixType::Zero(d.rows(), 1);
	i = 0;
	ct = 0;
	while (true) {
		d.row(ct) = directions.row(i);
		d.row(ct + 1) = -1 * directions.row(i);
		ex(ct) = extrema(idx.at(i));
		ex(ct + 1) = ex(ct);

		MatrixType::Index id;
		double tmp = (directions * d.row(ct + 1).transpose()).maxCoeff(&id);
		if (tmp > 0.95) {
			m.remove_row(directions, id);
			idx.erase(idx.begin() + id);
		}
		i += 1;
		ct += 2;

		if (i > directions.rows() - 1)
			break;
	}

	MatrixType ma = d;
	cout << " matrix rows,cols " << ma.rows() << " " << ma.cols() << endl << ma << endl;

	return 0;
}
