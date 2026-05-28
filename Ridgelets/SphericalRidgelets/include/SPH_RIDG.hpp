#include "SPH_RIDG.h"

#ifndef SPH_RIDG_IMPL
#define SPH_RIDG_IMPL

//Default Constructor
template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::SPH_RIDG()
{
	J = 2;
	rho = static_cast<pT>(1.0) / USER_RHO_DEFAULT;
	UM = UtilMath<pT, MT, VT>();
	init();
}

template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::~SPH_RIDG() {}

template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::SPH_RIDG(unsigned JIn, pT rhoIn) {
	J = JIn;
	rho = rhoIn;
	UM = UtilMath<pT, MT, VT>();
	init();
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::init() {
	// At the finest scale J the Gaussian/Weierstrass kernel is evaluated at
	// approximately l / 2^J, so its spectral tail behaves like
	// exp(-rho * l^2 / 4^J).  Truncate where that tail reaches 1e-6.
	mcut = (int)ceil(sqrt(-log(TAIL_THRESHOLD) * pow(4, J) / rho));
	mcut = mcut + mcut % 2;

	h = MT::Zero(mcut + 1, J + 1);
	h.col(0) = VT::LinSpaced(mcut + 1, 0, mcut);

	for (int i = 1; i < J + 1; ++i)
		h.col(i) = h.col(i - 1) * 0.5;

	h = (-rho * h.array().cwiseProduct(h.array() + 1.0)).exp();

	Lmd = MT::Zero(mcut + 1, 1);
	UM.fura(Lmd, mcut);
	psi = h;
	psi.col(0) = psi.col(0).cwiseProduct(Lmd);
	for (int i = 1; i < J + 1; ++i) {
		psi.col(i) = (h.col(i) - h.col(i - 1)).cwiseProduct(Lmd);
	}

	C = MT::Zero(mcut + 1, 1);
	C.col(0) = (2.0 * VT::LinSpaced(mcut + 1, 0, mcut).array() + 1.0) / (4 * UM.PI);

	t = ((psi.cwiseProduct(psi)).transpose() * C).array().sqrt();
	for (int i = 0; i < J + 1; ++i)
		psi.col(i) = psi.col(i) / t(i);

	// scale_weights(j) is the aligned ODF-domain response of one normalized
	// ridgelet atom at scale j.  RBasis uses sum_l C_l psi_lj P_l(mu); applying
	// the same Funk-Radon convention as QBasis evaluates P_l(0) at alignment.
	// UtilMath::fura stores the same P_l(0) sequence used by QBasis, up to the
	// repository's omitted global 2*pi factor, so these weights match -odf's
	// normalization convention.
	MT P0 = MT::Ones(1, mcut + 1);
	MT x0 = MT::Zero(1, 1);
	UM.polyleg(P0, x0, mcut);

	scale_weights = MT::Zero(J + 1, 1);
	for (int i = 0; i < J + 1; ++i) {
		MT r = C.cwiseProduct(psi.col(i));
		scale_weights(i, 0) = (P0 * r)(0, 0);
	}

	tau = 4.0 * log(10.0) / rho;
	m0 = (int)floor((-1.0 + sqrt(1.0 + 4.0 * tau)) / 2.0);
	M0 = VectorXi::Zero(J + 1, 1);

	for (int i = 0; i < J + 1; ++i)
		M0(i, 0) = (int)pow((pow(2, i) * m0 + 1), 2);
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::RBasis(MT& A, MT& u) {
	A.resize(u.rows(), M0.sum());
	A.setZero();

	MT P;
	MT X;
	MT x;

	MT v;
	MT vv;
	MT r;
	int I = 0;
	for (int i = 0; i < J + 1; ++i) {
		int K = M0(i);
		int N = 2 * K;
		v = MT::Zero(N, 3);
		UM.spiralsample(v, 2, N);
		vv = MT::Zero(K, 3);
		vv = v.topRows(K);

		r = C.cwiseProduct(psi.col(i));
		P = MT::Ones(u.rows(), mcut + 1);
		X = u * vv.transpose();
		for (int k = 0; k < K; ++k) {
			x = X.col(k);
			UM.polyleg(P, x, mcut);
			A.col(k + I) = P * r;
		}
		I += K;
	}
}

template<class pT, class MT, class VT>
const MT& SPH_RIDG<pT, MT, VT>::getScaleWeights() const {
	return scale_weights;
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::QBasis(MT& Q, MT& u) {
	Q.resize(u.rows(), M0.sum());
	Q.setZero();

	MT P;
	MT x;

	MT v;
	MT vv;
	MT r;
	MT CL = C.cwiseProduct(Lmd);

	int I = 0;
	for (int i = 0; i < J + 1; ++i) {
		int K = M0(i);
		int N = 2 * K;
		v = MT::Zero(N, 3);
		UM.spiralsample(v, 2, N);
		vv = v.topRows(K);

		r = CL.cwiseProduct(psi.col(i));
		for (int k = 0; k < K; ++k) {
			x = u * vv.row(k).transpose();
			unsigned Nr = (unsigned)x.rows();
			P = MT::Ones(Nr, mcut + 1);
			UM.polyleg(P, x, mcut);

			Q.col(k + I) = P * r;
		}
		I += K;
	}
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::FindMaxRidgeletMaxInDMRI(MT& fin, MT& R, pT nms_angle_degrees) {
	const int expected_coeffs = M0.sum();
	if (R.rows() != expected_coeffs)
		throw std::logic_error("Ridgelet coefficient count does not match the SPH_RIDG basis.");

	const bool use_abs_coefficients = false;
	const unsigned max_dirs = 6;
	const unsigned max_axes = max_dirs / 2;
	const pT nms_angle = std::max(static_cast<pT>(0.0), std::min(nms_angle_degrees, static_cast<pT>(90.0)));
	const pT nms_cos = std::cos(nms_angle * UM.PI / static_cast<pT>(180.0));

	struct Candidate {
		pT score;
		pT x;
		pT y;
		pT z;
	};

	vector<Candidate> candidates_template;
	candidates_template.reserve(expected_coeffs);

	for (int j = 0; j < J + 1; ++j) {
		int K = M0(j);
		int N = 2 * K;
		MT v = MT::Zero(N, 3);
		UM.spiralsample(v, 2, N);

		pT weight = scale_weights(j, 0);
		// The current FRT convention gives positive aligned responses.  Using
		// abs(weight) keeps the score invariant to an equivalent global sign
		// convention while the default coefficient rule stays positive-part.
		weight = std::abs(weight);

		for (int k = 0; k < K; ++k) {
			Candidate c;
			c.score = weight;
			c.x = v(k, 0);
			c.y = v(k, 1);
			c.z = v(k, 2);
			candidates_template.push_back(c);
		}
	}

	fin.resize((max_dirs * 3) + max_dirs, R.cols());
	fin.setZero((max_dirs * 3) + max_dirs, R.cols());

#if defined(_OPENMP)
	#pragma omp parallel for
#endif
	for (int i = 0; i < R.cols(); ++i)
	{
		vector<Candidate> candidates;
		candidates.reserve(candidates_template.size());

		int offset = 0;
		for (int j = 0; j < J + 1; ++j) {
			int K = M0(j);
			for (int k = 0; k < K; ++k) {
				pT coeff = R(offset + k, i);
				pT score;
				if (use_abs_coefficients)
					score = candidates_template[offset + k].score * std::abs(coeff);
				else
					score = candidates_template[offset + k].score * std::max(coeff, static_cast<pT>(0.0));

				if (score > static_cast<pT>(0.0)) {
					Candidate c = candidates_template[offset + k];
					c.score = score;
					candidates.push_back(c);
				}
			}
			offset += K;
		}

		std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
			return a.score > b.score;
		});

		vector<Candidate> selected;
		selected.reserve(max_axes);
		for (unsigned c = 0; c < candidates.size(); ++c) {
			bool suppressed = false;
			for (unsigned s = 0; s < selected.size(); ++s) {
				pT dot = candidates[c].x * selected[s].x +
					candidates[c].y * selected[s].y +
					candidates[c].z * selected[s].z;
				// Spherical ridgelet directions are antipodally symmetric, so
				// v and -v are the same fiber axis for non-maximum suppression.
				if (std::abs(dot) >= nms_cos) {
					suppressed = true;
					break;
				}
			}
			if (!suppressed) {
				selected.push_back(candidates[c]);
				if (selected.size() == max_axes)
					break;
			}
		}

		for (unsigned a = 0; a < selected.size(); ++a) {
			unsigned row = 8 * a;
			fin(row, i) = selected[a].x;
			fin(row + 1, i) = selected[a].y;
			fin(row + 2, i) = selected[a].z;
			fin(row + 3, i) = selected[a].score;
			fin(row + 4, i) = -selected[a].x;
			fin(row + 5, i) = -selected[a].y;
			fin(row + 6, i) = -selected[a].z;
			fin(row + 7, i) = selected[a].score;
		}
	}
}

template<class pT, class MT, class VT>
bool SPH_RIDG<pT, MT, VT>::TestDirectRidgeletMaxima(pT& abs_dot, pT& score, int& coeff_index) {
	int offset = 0;
	for (int j = 0; j < J; ++j)
		offset += M0(j);

	coeff_index = offset;
	MT R = MT::Zero(M0.sum(), 1);
	R(coeff_index, 0) = static_cast<pT>(1.0);

	MT fin;
	FindMaxRidgeletMaxInDMRI(fin, R);

	int K = M0(J);
	int N = 2 * K;
	MT v = MT::Zero(N, 3);
	UM.spiralsample(v, 2, N);

	abs_dot = std::abs(fin(0, 0) * v(0, 0) +
		fin(1, 0) * v(0, 1) +
		fin(2, 0) * v(0, 2));
	score = fin(3, 0);

	return abs_dot > static_cast<pT>(1.0 - 1e-6) && score > static_cast<pT>(0.0);
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::normBasis(MT& mat) {
	MT e = mat * mat.transpose();
	SelfAdjointEigenSolver<MT> eigensolver(mat.rows());
	eigensolver.compute(e, EigenvaluesOnly);
	pT lVal = eigensolver.eigenvalues()[mat.rows() - 1];
	mat = (1 / sqrt(lVal)) * mat;
}

#endif
