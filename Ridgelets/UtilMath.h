#ifndef UTILMATH_H
#define UTILMATH_H

#include <iostream>
#include "rdgls_types.h"

using namespace std;

class UtilMath {
public:
	//Useful constants
	const double PI = atan(1.0) * 4.0;

	//Functions
	void spiralsample(MatrixType& u, unsigned flg, unsigned N);
	void fura(MatrixType& Lmd, unsigned n);
	void polyleg(MatrixType& P, MatrixType x, unsigned n);
	MatrixType convhulln(MatrixType & u);
	void unique_rows(vector<int>& uniques, MatrixType & U);
	void unique_sorted(vector<unsigned>& uniques, MatrixType & U);
	void ind_sort(MatrixType & matrix, multimap<double, unsigned>& indx, unsigned col_n);
	void column_find(std::vector<Eigen::Index>& index, MatrixType & arr, unsigned col_n, bool equal, int val);
	void icosahedron(MatrixType& u,MatrixType & faces, unsigned level);
	void index_and_flat(MatrixType & u, vector<Eigen::Index>& a, MatrixType & fcs, unsigned sz, unsigned col);
	void FindConnectivity(vector<vector<unsigned>>& conn, MatrixType & fcs, unsigned N);
};

#endif // ! UTILMATH_H