#ifndef UTILMATH_H
#define UTILMATH_H

#include "rdgls_types.h"

using namespace std;

class UtilMath {
public:
	struct cmp_v
	{
    		const vector<precisionType> & value_vector;

    		cmp_v(const vector<precisionType> & val_vec):
        		value_vector(val_vec) {}

    		bool operator()(size_t i1, size_t i2)
    		{   
        		return value_vector[i1] < value_vector[i2];
    		}
	};

	// Useful constants
	static const precisionType PI;

	UtilMath();
	~UtilMath();

	// Functions
	void spiralsample(MatrixType & u, unsigned flg, unsigned N);
	void fura(MatrixType & Lmd, unsigned n);
	void polyleg(MatrixType & P, MatrixType & x, unsigned n);
	MatrixType convhull3_1(MatrixType & u);
	void unique_rows(vector<int> & uniques, MatrixType & U);
	void unique_sorted(vector<unsigned> & uniques, MatrixType & U);
	void ind_sort(MatrixType & matrix, multimap<precisionType, unsigned> & indx, unsigned col_n);
	void ind_sort_vec(MatrixType & vec, vector<size_t> & indx);
	void column_find(std::vector<Eigen::Index>& index, MatrixType & arr, unsigned col_n, bool equal, int val);
	void icosahedron(MatrixType& u,MatrixType & faces, unsigned level);
	void index_and_flat(MatrixType & u, vector<Eigen::Index>& a, MatrixType & fcs, unsigned sz, unsigned col);
	void FindConnectivity(vector<vector<unsigned>>& conn, MatrixType & fcs, unsigned N);
	void remove_row(MatrixType & a, MatrixType::Index del);
	void FindODFMaxima(MatrixType & ex, MatrixType & d, MatrixType & W, vector<vector<unsigned>>& conn, MatrixType & u, float thresh);
	void FindMaxODFMaxInDMRI(MatrixType & fin, MatrixType & Q, MatrixType & C, vector<vector<unsigned>>& conn, MatrixType & nu, float thresh);
};

#endif // ! UTILMATH_H
