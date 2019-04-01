#include "UtilMath.h"

#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"

UtilMath::UtilMath() {}
UtilMath::~UtilMath() { cout << "UtilMath desctructed" << endl; }

void UtilMath::spiralsample(MatrixType& u, unsigned flg, unsigned N)
{
	/*
	Spiral sampling on the sphere
	*/
	//z=(1-1/N:-2/N:1/N-1)';
	MatrixType z(N, 1);
	double val = 1.0 - (1.0 / N);
	for (unsigned int i = 0; i < N; ++i, val += (-2.0 / N))
		z(i) = val;

	//r=sqrt(1-z.*z);
	MatrixType r(N, 1);
	r = (1.0 - z.cwiseProduct(z).array()).sqrt();

	MatrixType Long(N, 1);
	switch (flg)
	{
	case 1:
	{
		// Long=[0; cumsum((3.6/sqrt(N))./r(1:N-1))];
		double sqrtN = sqrt(N);
		Long(0) = 0.0;
		for (unsigned i = 1; i < N; ++i)
			Long(i) = Long(i - 1) + ((3.6 / sqrtN) / r(i - 1));
	}
	break;
	case 2:
	{
		// Long=(pi*(3-sqrt(5)))*(0:N-1)';
		double pref_const = UtilMath::PI * (3.0 - sqrt(5.0));
		Long = pref_const * VectorXd::LinSpaced(N, 0, N - 1);
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

void UtilMath::fura(MatrixType& Lmd, unsigned n)
{
	//from fura.m

	//Lmd=ones(n+1,1);

	//Lmd(2:2:n+1)=0;
	Lmd(0) = 1.0;
	double c;

	//for k=2:2:n
	for (unsigned k = 2; k < n + 1; k += 2) {
		//Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
		c = (double)k;
		Lmd.row(k) = -((c - 1.0) / c)*Lmd.row(k - 2);
	}

	//Lmd=(2*pi)*Lmd;
	//Lmd = Lmd * (2 * PI);
}

void UtilMath::polyleg(MatrixType& P, MatrixType& x, unsigned n)
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
			double c1 = (2.0 * k - 1.0) / k;
			double c2 = (k - 1.0) / k;
			P.col(k) = c1 * x.cwiseProduct(P.col(k - 1)) - c2 * P.col(k - 2);
		}
	}
}

MatrixType UtilMath::convhull3_1(MatrixType& u) {
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

	MatrixType faces(nFaces, 3);
	for (int i = 0; i < nFaces; i++)
		for (int j = 0; j < 3; j++)
			faces(i, j) = *faceIndices++;

	free(vertices);

	return faces;
}

void UtilMath::unique_rows(vector<int>& uniques, MatrixType& U) {
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
	size_t added_length = nCols * to_string(U(0, 0)).length();
	key.reserve(key.length() + added_length);

	// Iterate over matrix
	for (unsigned i = 0; i < U.rows(); ++i) {
		// Create unique key from row valuse
		key.clear();
		for (unsigned j = 0; j < nCols; ++j) {
			key.append(to_string(U(i, j)));
		}

		// If element not exists in hash table
		if (hTable.count(key) == 0) {
			hTable.insert(pair<string, bool>(key, true));
			uniques.push_back(i);
		}
	}
}

void UtilMath::unique_sorted(vector<unsigned>& uniques, MatrixType& U) {
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

void UtilMath::ind_sort(MatrixType& matrix, multimap<double, unsigned>& indx, unsigned col_n) {
	/*
	Return indexies indx in ascending order of sorted column col_n of the matrix
	*/
	// Matrix column to std vector
	vector<double> uc3;
	unsigned orig_size = matrix.col(col_n).size();
	uc3.resize(orig_size);
	VectorXd::Map(&uc3[0], orig_size) = matrix.col(col_n);

	// Mapping from value to index and so make sort ascending
	for (auto it = uc3.begin(); it != uc3.end(); ++it)
		indx.insert(make_pair(*it, it - uc3.begin()));
}

void UtilMath::ind_sort_vec(MatrixType& vec, vector<size_t> & indx) {
	/*
	Return indexies indx in descending order of sorted vector
	*/
	// Matrix to std vector
	vector<double> uc3;
	unsigned orig_size = vec.size();
	uc3.resize(orig_size);
	VectorXd::Map(&uc3[0], orig_size) = vec;

	// initialize original index locations
	indx.resize(uc3.size());
	iota(indx.begin(), indx.end(), 0);

	// sort indexes based on comparing values in v
	sort(indx.rbegin(), indx.rend(), [&uc3](size_t i1, size_t i2) {return uc3[i1] < uc3[i2]; });
}

void UtilMath::column_find(std::vector<Eigen::Index>& index, MatrixType& arr, unsigned col_n, bool equal, int val) {
	/*
	Looking for columns of MatrixType matrix. bool equal is basically to comply with matlab notation of
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

void UtilMath::icosahedron(MatrixType& u, MatrixType& faces, unsigned level) {
	/*
	Build icosahedron based sampling array and their respective faces
	*/
	cout << "Start computing icosahedron..." << endl;
	double C = 1 / sqrt(1.25);
	MatrixType t = (2 * PI / 5.0) * VectorXd::LinSpaced(5, 0, 4);
	MatrixType u1(5, 3);
	u1 << C * t.array().cos(), C * t.array().sin(), C * 0.5 * MatrixType::Ones(5, 1);
	MatrixType u2(5, 3);
	u2 << C * (t.array() + 0.2 * PI).cos(), C * (t.array() + 0.2 * PI).sin(), -0.5 * C * MatrixType::Ones(5, 1);
	u.resize(12, 3);
	u << 0, 0, 1, u1, u2, 0, 0, -1;
	MatrixType u_final;

	if (level > 0) {
		for (unsigned lev = 1; lev <= level; ++lev) {
			MatrixType fcs = convhull3_1(u);

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
			unique_rows(uniques, U);

			// Normalize and add to u
			unsigned u_len = u.rows();
			unsigned unique_len = uniques.size();
			u.conservativeResize(u.rows() + unique_len, u.cols());

			for (unsigned i = 0, j = 0 + u_len; i < unique_len, j < unique_len + u_len; ++i, ++j)
				u.row(j) = U.row(uniques.at(i)) / U.row(uniques.at(i)).norm();
		}

		// Sorting u by 3rd col
		multimap<double, unsigned> ind;
		ind_sort(u, ind, 2);

		// Using indicies in reverse order gives us desired descending order
		unsigned i = 0;
		MatrixType u_sorted(u.rows(), 3);
		for (multimap<double, unsigned>::reverse_iterator it = ind.rbegin(); it != ind.rend(); ++it) {
			u_sorted.row(i) = u.row(it->second);
			++i;
		}

		// Find indicies where 3rd column eq 0
		std::vector<Eigen::Index> index;
		column_find(index, u_sorted, 2, true, 0);

		// v matrix part of u where 3rd col eq 0
		unsigned N_index = index.size();
		MatrixType v(N_index, 3);
		for (unsigned i = 0; i < N_index; ++i)
			v.row(i) = u_sorted.row(index.at(i));

		// Sort v by 2nd column
		multimap<double, unsigned> ind_v;
		ind_sort(v, ind_v, 1);

		// Using indicies in reverse order gives us desired descending order
		i = 0;
		for (multimap<double, unsigned>::reverse_iterator it = ind_v.rbegin(); it != ind_v.rend(); ++it) {
			u_sorted.row(index.at(i)) = v.row(it->second);
			++i;
		}
		u = u_sorted;
	}

	// Normalize
	u = u.array().colwise() / (u.rowwise().norm().array() + 2.2204e-16);

	faces = convhull3_1(u);
}

void UtilMath::index_and_flat(MatrixType& u, vector<Eigen::Index>& a, MatrixType& fcs, unsigned sz, unsigned col) {
	/*
	Create sub array of fcs of size sz starting in column col and rows a
	*/
	for (unsigned i = 0; i < a.size(); ++i)
		u.row(i) = fcs.block(a.at(i), col, 1, sz);
	u.resize(sz * a.size(), 1);
}

void UtilMath::FindConnectivity(vector<vector<unsigned>>& conn, MatrixType& fcs, unsigned N) {
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

		MatrixType u1(a1.size(), 2);
		index_and_flat(u1, a1, fcs, 2, 1);

		MatrixType u2(a2.size(), 1);
		index_and_flat(u2, a2, fcs, 1, 0);

		MatrixType u3(a2.size(), 1);
		index_and_flat(u3, a2, fcs, 1, 2);

		MatrixType u4(a3.size(), 2);
		index_and_flat(u4, a3, fcs, 2, 0);

		MatrixType u_out(u1.size() + u2.size() + u3.size() + u4.size(), 1);
		u_out << u1, u2, u3, u4;

		vector<unsigned> un;
		unique_sorted(un, u_out);

		conn.push_back(un);
	}
}

void UtilMath::remove_row(MatrixType& a, MatrixType::Index del)
{
	unsigned cols = a.cols();
	unsigned rows = a.rows() - 1;

	if (del < rows)
		a.block(del, 0, rows - del, cols) = a.block(del + 1, 0, rows - del, cols);

	a.conservativeResize(rows, cols);
}

void UtilMath::FindODFMaxima(MatrixType& ex, MatrixType& d, MatrixType& W,
	vector<vector<unsigned>>& conn, MatrixType& u, float thresh)
{
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
					double maxw = 0;
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
	MatrixType directions(u_extrema_length, 3);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		directions.row(i) = u.row(u_extrema.at(i));

	//sort(W(extrema),'descend')
	MatrixType W_e(u_extrema_length, 1);
	for (unsigned i = 0; i < u_extrema_length; ++i)
		W_e(i) = W(u_extrema.at(i));

	vector<size_t> idx;
	ind_sort_vec(W_e, idx);

	MatrixType directions_sorted(directions.rows(), directions.cols());
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

		MatrixType::Index id;
		double tmp = (directions_sorted * d.row(ct).transpose()).maxCoeff(&id);
		if (tmp > 0.95) {
			remove_row(directions_sorted, id);
			idx.erase(idx.begin() + id);
		}
		i += 1;
		ct += 2;

		if (i > directions_sorted.rows() - 1)
			break;
	}
}

void UtilMath::FindMaxODFMaxInDMRI(MatrixType& fin, MatrixType& cnt, MatrixType& Q, 
	MatrixType& C, vector<vector<unsigned>>& conn, MatrixType& nu, float thresh)
{
	fin.resize((6 * 3) + 6, C.cols());
	fin.setZero((6 * 3) + 6, C.cols());

	cnt.resize(1, C.cols());
	cnt.setZero(1, C.cols());

	#pragma omp parallel for
	for (int i = 0; i < C.cols(); ++i)
	{
		MatrixType exe_vol;
		MatrixType dir_vol;
		MatrixType vol = Q * C.col(i);

		FindODFMaxima(exe_vol, dir_vol, vol, conn, nu, thresh);

		if (exe_vol.rows() < 16)
			cnt(i) = exe_vol.rows();
		else
			cnt(i) = 0;

		unsigned ex_sz = std::min((int)exe_vol.rows(), 6);
		dir_vol.conservativeResize(dir_vol.rows() * 3, 1);

		// Loop over to make array in a form (x y z val) for direction (max 6)
		for (unsigned j = 0, k = 0, z = 0; j < ex_sz; ++j, k += 3, z += 4)
		{
			fin(z, i) = dir_vol(k + 1);
			fin(z + 1, i) = dir_vol(k + 2);
			fin(z + 2, i) = dir_vol(k + 3);
			fin(z + 3, i) = vol(exe_vol(j));
		}
	}
}