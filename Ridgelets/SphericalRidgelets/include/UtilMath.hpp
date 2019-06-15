#include "UtilMath.h"

#ifndef UTILMATH_IMPL
#define UTILMATH_IMPL

template <class pT, class MT, class VT>
const pT UtilMath<pT, MT, VT>::PI = std::atan(1.0) * 4.0;

template <class pT, class MT, class VT>
UtilMath<pT, MT, VT>::UtilMath() {}

template <class pT, class MT, class VT>
UtilMath<pT, MT, VT>::~UtilMath() {}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::spiralsample(MT& u, unsigned flg, unsigned N)
{
	/*
	Spiral sampling on the sphere
	*/
	//z=(1-1/N:-2/N:1/N-1)';
	MT z(N, 1);
	pT val = 1.0 - (1.0 / N);
	for (unsigned int i = 0; i < N; ++i, val += (-2.0 / N))
		z(i) = val;

	//r=sqrt(1-z.*z);
	MT r(N, 1);
	r = (1.0 - z.cwiseProduct(z).array()).sqrt();

	MT Long(N, 1);
	switch (flg)
	{
	case 1:
	{
		// Long=[0; cumsum((3.6/sqrt(N))./r(1:N-1))];
		pT sqrtN = sqrt(N);
		Long(0) = 0.0;
		for (unsigned i = 1; i < N; ++i)
			Long(i) = Long(i - 1) + ((3.6 / sqrtN) / r(i - 1));
	}
	break;
	case 2:
	{
		// Long=(pi*(3-sqrt(5)))*(0:N-1)';
		pT pref_const = UtilMath::PI * (3.0 - sqrt(5.0));
		Long = pref_const * VT::LinSpaced(N, 0, N - 1);
	}
	break;
	default:
		cerr << "Invalid Sampling Option \n";
		throw;
	}

	//u=[r.*cos(long) r.*sin(long) z];
	//u=u./repmat(sqrt(sum(u.^2,2)),[1 3]);

	u.col(0) = r.array() * Long.array().cos();
	u.col(1) = r.array() * Long.array().sin();
	u.col(2) = z.array();
	for (unsigned i = 0; i < N; ++i)
		u.row(i).normalize();
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::fura(MT& Lmd, unsigned n)
{
	//from fura.m

	//Lmd=ones(n+1,1);

	//Lmd(2:2:n+1)=0;
	Lmd(0) = 1.0;
	pT c;

	//for k=2:2:n
	for (unsigned k = 2; k < n + 1; k += 2) {
		//Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
		c = (pT)k;
		Lmd.row(k) = -((c - 1.0) / c)*Lmd.row(k - 2);
	}

	//Lmd=(2*pi)*Lmd;
	//Lmd = Lmd * (2 * PI);
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::polyleg(MT& P, MT& x, unsigned n)
{
	//from polyleg.m

	if (x.cwiseAbs().maxCoeff() > 1) {
		cerr << "Values should be in range [-1, 1] \n";
		throw;
	}

	switch (n)
	{
	case 0:
		break;
	case 1:
		P.col(1) = x;
		break;
	default:
		P.col(1) = x;
		for (unsigned k = 2; k < n + 1; ++k)
		{
			pT c1 = (2.0 * k - 1.0) / k;
			pT c2 = (k - 1.0) / k;
			P.col(k) = c1 * x.cwiseProduct(P.col(k - 1)) - c2 * P.col(k - 2);
		}
	}
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::convhull3_1(MT& u, MT& fcs) {
	unsigned n = u.rows();
	ch_vertex* vertices;
	vertices = (ch_vertex*)malloc(n * sizeof(ch_vertex));
	for (unsigned i = 0; i < n; ++i) {
		vertices[i].x = u(i, 0);
		vertices[i].y = u(i, 1);
		vertices[i].z = u(i, 2);
	}

	int* faceIndices = NULL;
	int nFaces;
	convhull_3d_build(vertices, n, &faceIndices, &nFaces);

	fcs.resize(nFaces, 3);
	for (int i = 0; i < nFaces; i++)
		for (int j = 0; j < 3; j++)
			fcs(i, j) = *faceIndices++;

	free(vertices);
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::unique_rows(vector<int>& uniques, MT& U) {
	/*
		Find unique rows
		Not the fanciest and most optimal way, but
		uses kind of hashtable as serious boys usually do :)
	*/

	unsigned nCols = U.cols();

	// Define hashtable
	unordered_map<string, bool> hTable;

	// Preallocate string for faster string concatenation
	string key;
	// long double always for to_string function
	size_t added_length = nCols * to_string(static_cast<long double>(U(0, 0))).length();

	key.reserve(key.length() + added_length);

	// Iterate over matrix
	for (unsigned i = 0; i < U.rows(); ++i) {
		// Create unique key from row values
		key.clear();
		for (unsigned j = 0; j < nCols; ++j) {
			key.append(to_string(static_cast<long double>(U(i, j))));
		}

		// If element not exists in hash table
		if (hTable.count(key) == 0) {
			hTable.insert(pair<string, bool>(key, true));
			uniques.push_back(i);
		}
	}
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::unique_sorted(vector<unsigned>& uniques, MT& U) {
	/*
		Find unique values and return them sorted (same as default matlab version)
		Now working with 1D and integer values of U only
	*/
	// Define hashtable
	unordered_map<long, bool> hTable;

	// Iterate over matrix
	for (unsigned i = 0; i < U.size(); ++i) {
		// If element not exists in hash table
		if (hTable.count(U(i)) == 0) {
			hTable.insert(pair<long, bool>(U(i), true));
			uniques.push_back(U(i));
		}
	}
	sort(uniques.begin(), uniques.end());
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::ind_sort(MT& matrix, multimap<pT, unsigned>& indx, unsigned col_n) {
	/*
	Return indexies indx in ascending order of sorted column col_n of the matrix
	*/
	// Matrix column to std vector
	vector<pT> uc3;
	unsigned orig_size = matrix.col(col_n).size();
	uc3.resize(orig_size);
	VT::Map(&uc3[0], orig_size) = matrix.col(col_n);

	// Mapping from value to index and so make sort ascending
	for (auto it = uc3.begin(); it != uc3.end(); ++it)
		indx.insert(make_pair(*it, it - uc3.begin()));
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::ind_sort_vec(MT& vec, vector<size_t> & indx) {
	/*
	Return indexies indx in descending order of sorted vector
	*/
	// Matrix to std vector
	vector<pT> uc3;
	unsigned orig_size = vec.size();
	uc3.resize(orig_size);
	VT::Map(&uc3[0], orig_size) = vec;

	// initialize original index locations
	indx.resize(uc3.size());
	iota(indx.begin(), indx.end(), 0);
	sort(indx.rbegin(), indx.rend(), cmp_v(uc3));
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::column_find(std::vector<Eigen::Index>& index, MT& arr, unsigned col_n, bool equal, int val) {
	/*
	Looking for columns of MT matrix. bool equal is basically to comply with matlab notation of
	== if true and ~= if false. So you can check if col_n equal or not to val
	*/
	for (Eigen::Index i = 0; i < arr.rows(); ++i) {
		if (equal) {
			if (arr.col(col_n)(i) == val)
				index.push_back(i);
		}
		else {
			if (arr.col(col_n)(i) != val)
				index.push_back(i);
		}
	}
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::icosahedron(MT& u, MT& faces, unsigned level) {
	/*
	Build icosahedron based sampling array and their respective faces
	*/
	cout << "Start computing icosahedron..." << endl;
	pT C_init = 1 / sqrt(1.25);
	MT t = (2 * PI / 5.0) * VT::LinSpaced(5, 0, 4);
	MT u1(5, 3);
	u1 << C_init * t.array().cos(), C_init * t.array().sin(), C_init * 0.5 * MT::Ones(5, 1);
	MT u2(5, 3);
	u2 << C_init * (t.array() + 0.2 * PI).cos(), C_init * (t.array() + 0.2 * PI).sin(), -0.5 * C_init * MT::Ones(5, 1);
	u.resize(12, 3);
	u << 0, 0, 1, u1, u2, 0, 0, -1;

	MT U;
	MT A;
	MT B;
	MT C;
	MT fcs;
	vector<int> uniques;

	if (level > 0) {
		for (unsigned lev = 1; lev <= level; ++lev) {
			convhull3_1(u, fcs);

			unsigned N = fcs.rows();
			U = MT::Zero(3 * N, 3);
			for (unsigned k = 0; k < N; ++k) {
				A = u.row(fcs(k, 0));
				B = u.row(fcs(k, 1));
				C = u.row(fcs(k, 2));
				U.template block<3, 3>(3 * k, 0) << 0.5 * (A + B), 0.5 * (B + C), 0.5 * (A + C);
			}

			unique_rows(uniques, U);

			// Normalize and add to u
			unsigned u_len = u.rows();
			unsigned unique_len = uniques.size();

			u.conservativeResize(u.rows() + unique_len, u.cols());

			for (unsigned i = 0, j = 0 + u_len; j < unique_len + u_len; ++i, ++j)
				u.row(j) = U.row(uniques.at(i)) / U.row(uniques.at(i)).norm();
		}

		// Sorting u by 3rd col
		multimap<pT, unsigned> ind;
		ind_sort(u, ind, 2);

		// Using indicies in reverse order gives us desired descending order
		unsigned i = 0;
		MT u_sorted(u.rows(), 3);
		for (typename multimap<pT, unsigned>::reverse_iterator it = ind.rbegin(); it != ind.rend(); ++it) {
			u_sorted.row(i) = u.row(it->second);
			++i;
		}

		// Find indicies where 3rd column eq 0
		std::vector<Eigen::Index> index;
		column_find(index, u_sorted, 2, true, 0);

		// v matrix part of u where 3rd col eq 0
		unsigned N_index = index.size();
		MT v(N_index, 3);
		for (unsigned k = 0; k < N_index; ++k)
			v.row(k) = u_sorted.row(index.at(k));

		// Sort v by 2nd column
		multimap<pT, unsigned> ind_v;
		ind_sort(v, ind_v, 1);

		// Using indicies in reverse order gives us desired descending order
		i = 0;
		for (typename multimap<pT, unsigned>::reverse_iterator it = ind_v.rbegin(); it != ind_v.rend(); ++it) {
			u_sorted.row(index.at(i)) = v.row(it->second);
			++i;
		}
		u = u_sorted;
	}

	// Normalize
	u = u.array().colwise() / (u.rowwise().norm().array() + 2.2204e-16);

	convhull3_1(u, faces);
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::index_and_flat(MT& u, vector<Eigen::Index>& a, MT& fcs, unsigned sz, unsigned col) {
	/*
	Create sub array of fcs of size sz starting in column col and rows a
	*/
	for (unsigned i = 0; i < a.size(); ++i)
		u.row(i) = fcs.block(a.at(i), col, 1, sz);
	u.resize(sz * a.size(), 1);
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::FindConnectivity(vector<vector<unsigned>>& conn, MT& fcs, unsigned N) {
	/*
	FindConnectivity function. It is return vector of vectors because in matlab function it returns dynamic size
	cell array. So it is the best way I found to store arrays of indcicies with different sizes.
	*/
	for (unsigned i = 0; i < N; ++i) {
		vector<Eigen::Index> a1;
		column_find(a1, fcs, 0, true, i);

		vector<Eigen::Index> a2;
		column_find(a2, fcs, 1, true, i);

		vector<Eigen::Index> a3;
		column_find(a3, fcs, 2, true, i);

		MT u1(a1.size(), 2);
		index_and_flat(u1, a1, fcs, 2, 1);

		MT u2(a2.size(), 1);
		index_and_flat(u2, a2, fcs, 1, 0);

		MT u3(a2.size(), 1);
		index_and_flat(u3, a2, fcs, 1, 2);

		MT u4(a3.size(), 2);
		index_and_flat(u4, a3, fcs, 2, 0);

		MT u_out(u1.size() + u2.size() + u3.size() + u4.size(), 1);
		u_out << u1, u2, u3, u4;

		vector<unsigned> un;
		unique_sorted(un, u_out);

		conn.push_back(un);
	}
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::remove_row(MT& a, Eigen::Index del)
{
	unsigned cols = a.cols();
	unsigned rows = a.rows() - 1;

	if (del < rows)
		a.block(del, 0, rows - del, cols) = a.block(del + 1, 0, rows - del, cols);

	a.conservativeResize(rows, cols);
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::FindODFMaxima(MT& ex, MT& d, VT& W,
	vector<vector<unsigned>>& conn, MT& u, pT thresh, unsigned& n_of_dirs)
{
	// Standart min-max normalization
	pT W_min = W.minCoeff();
	W = (W.array() - W_min) / (W.maxCoeff() - W_min);

	// Find maxima above this point
	MT used = MT::Zero(W.rows(), W.cols());

	// used(W <= thresh) = 1
	for (unsigned i = 0; i < W.size(); ++i)
		if (W(i) <= thresh)
			used(i) = 1;

	unsigned ct = 0;

	MT extrema(0, 0);

	for (unsigned n = 0; n < used.size(); ++n) {
		if (used(n) == 0) {
			unsigned j = n;
			bool reached_maxima = false;
			while (!reached_maxima) {
				// if (any(W(conn(j).elem) >= W(j)))
				bool if_any = false;
				unsigned conn_row_length = conn[j].size();
				for (unsigned i = 0; i < conn_row_length; i++) {
					if (W(conn[j][i]) >= W(j)) {
						if_any = true;
						break; // trick to speed up computations.
					}
				}
				if (if_any) {
					// [maxw id] = max(W(conn(j).elem))
					unsigned id = 0;
					pT maxw = 0;
					for (unsigned i = 0; i < conn_row_length; i++) {
						if (W(conn[j][i]) > maxw) {
							maxw = W(conn[j][i]);
							id = i;
						}
					}

					// We have already traveled this path
					if (used(conn[j][id]))
						reached_maxima = true;

					// used(conn(j).elem) = 1
					for (unsigned i = 0; i < conn_row_length; ++i) {
						used(conn[j][i]) = 1;
					}

					// j = conn(j).elem(id)
					j = conn[j][id];
				}
				else {
					reached_maxima = true;
					extrema.conservativeResize(1, extrema.cols() + 1);
					extrema(ct) = j;
					for (unsigned i = 0; i < conn_row_length; ++i) {
						used(conn[j][i]) = 1;
					}

					ct += 1;
				}
			}
		}
	}

	if (extrema.size() == 0) {
		extrema.conservativeResize(1, 1);
		extrema(0) = 1;
	}

	vector<unsigned> u_extrema;
	unique_sorted(u_extrema, extrema);

	unsigned u_extrema_length = u_extrema.size();
	MT directions(u_extrema_length, 3);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		directions.row(i) = u.row(u_extrema.at(i));

	//sort(W(extrema),'descend')
	MT W_e(u_extrema_length, 1);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		W_e(i) = W(u_extrema.at(i));

	vector<size_t> idx;
	ind_sort_vec(W_e, idx);

	MT directions_sorted(directions.rows(), directions.cols());
	for (unsigned i = 0; i < directions.rows(); ++i) {
		directions_sorted.row(i) = directions.row(idx.at(i));
	}

	unsigned extrema_size = extrema.size();
	if (extrema_size < 2) {
		extrema_size = 2;
	}

	d.resize((unsigned)ceil(extrema_size / 2.0) * 2, 3);
	ex.resize(d.rows(), 1);

	d.setZero((unsigned)ceil(extrema_size / 2.0) * 2, 3);
	ex.setZero(d.rows(), 1);

	unsigned i = 0;
	ct = 1;
	while (true) {
		if (ct > d.rows()) {
			d.conservativeResize(ct + 1, d.cols());
		}
		if (ct > ex.rows()) {
			ex.conservativeResize(ct + 1, 1);
		}
		d.row(ct - 1) = directions_sorted.row(i);
		d.row(ct) = -1 * directions_sorted.row(i);

		ex(ct - 1) = u_extrema.at(idx.at(i));
		ex(ct) = ex(ct - 1);

		typename MT::Index id;
		pT tmp = (directions_sorted * d.row(ct).transpose()).maxCoeff(&id);
		if (tmp > 0.95) {
			remove_row(directions_sorted, id);
			idx.erase(idx.begin() + id);
		}
		i += 1;
		ct += 2;

		if (i > directions_sorted.rows() - 1)
			break;
	}

	n_of_dirs = static_cast<unsigned>(ex.rows()) / 2;
}

template <class pT, class MT, class VT>
void UtilMath<pT, MT, VT>::FindMaxODFMaxInDMRI(MT& fin, MT& Q, MT& C,
	vector<vector<unsigned>>& conn, MT& nu, pT thresh)
{
	fin.resize((6 * 3) + 6, C.cols());
	fin.setZero((6 * 3) + 6, C.cols());

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (int i = 0; i < C.cols(); ++i)
	{
		MT exe_vol;
		MT dir_vol;
		VT vol = Q * C.col(i);
		VT dirs_vec;
		unsigned ND;

		FindODFMaxima(exe_vol, dir_vol, vol, conn, nu, thresh, ND);

		unsigned ex_sz = std::min((int)exe_vol.rows(), 6);

		dirs_vec.resize(dir_vol.rows() * 3);
		for (int y = 0, r = 0; r < dir_vol.rows(); y += 3, ++r) {
			dirs_vec(y) = dir_vol(r, 0);
			dirs_vec(y + 1) = dir_vol(r, 1);
			dirs_vec(y + 2) = dir_vol(r, 2);
		}

		// Loop over to make array in a form (x y z val) for direction (max 6)
		for (unsigned j = 0, k = 0, z = 0; j < ex_sz; ++j, k += 3, z += 4)
		{
			fin(z, i) = dirs_vec(k);
			fin(z + 1, i) = dirs_vec(k + 1);
			fin(z + 2, i) = dirs_vec(k + 2);
			fin(z + 3, i) = vol(exe_vol(j));
		}
	}
}

#endif
