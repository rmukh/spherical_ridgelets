#ifndef UTILMATH_H
#define UTILMATH_H

#include "main_externals.h"
#include "convhull_3d.h"

template <class pT, class MT, class VT>
class UtilMath {
public:
	struct cmp_v
	{
    		const vector<pT> & value_vector;

    		cmp_v(const vector<pT> & val_vec):
        		value_vector(val_vec) {}

    		bool operator()(size_t i1, size_t i2)
    		{   
        		return value_vector[i1] < value_vector[i2];
    		}
	};

	// Useful constants
	static const pT PI;

	UtilMath();
	~UtilMath();

	// Functions
	void spiralsample(MT & u, unsigned flg, unsigned N);
	void fura(MT & Lmd, unsigned n);
	void polyleg(MT & P, MT & x, unsigned n);
	void convhull3_1(MT & u, MT & fcs);
	void unique_rows(vector<int> & uniques, MT & U);
	void unique_sorted(vector<unsigned> & uniques, MT & U);
	void ind_sort(MT & matrix, multimap<pT, unsigned> & indx, unsigned col_n);
	void ind_sort_vec(MT & vec, vector<size_t> & indx);
	void column_find(std::vector<Eigen::Index>& index, MT & arr, unsigned col_n, bool equal, int val);
	void icosahedron(MT& u,MT & faces, unsigned level);
	void index_and_flat(MT & u, vector<Eigen::Index>& a, MT & fcs, unsigned sz, unsigned col);
	void FindConnectivity(vector<vector<unsigned>>& conn, MT & fcs, unsigned N);
	void remove_row(MT & a, Eigen::Index del);
	void FindODFMaxima(MT & ex, MT & d, VT & W, vector<vector<unsigned>>& conn, MT & u, pT thresh, unsigned & n_of_dirs);
	void FindMaxODFMaxInDMRI(MT & fin, MT & Q, MT & C, vector<vector<unsigned>>& conn, MT & nu, pT thresh);
};

#include "UtilMath.hpp"

#endif // ! UTILMATH_H
