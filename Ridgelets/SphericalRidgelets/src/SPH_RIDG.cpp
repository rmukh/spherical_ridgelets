#include "SPH_RIDG.h"

#ifndef SPH_RIDG_IMPL
#define SPH_RIDG_IMPL

//Default Constructor
template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::SPH_RIDG()
{
	J = 2;
	rho = 0.5;
	UM = UtilMath<pT, MT>();
	init();
}

template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::~SPH_RIDG() {}

template<class pT, class MT, class VT>
SPH_RIDG<pT, MT, VT>::SPH_RIDG(unsigned JIn, pT rhoIn) {
	J = JIn;
	rho = rhoIn;
	UM = UtilMath<pT, MT>();
	init();
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::init() {
	mcut = (int)ceil(sqrt(-log(1e-6) * pow(4, J) / rho));
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

	tau = 4.0 * log(10.0) / rho;
	m0 = (int)floor((-1.0 + sqrt(1.0 + 4.0 * tau)) / 2.0);
	M0 = VectorXi::Zero(J + 1, 1);

	for (int i = 0; i < J + 1; ++i)
		M0(i, 0) = (int)pow((pow(2, i) * m0 + 1), 2);
}

template<class pT, class MT, class VT>
void SPH_RIDG<pT, MT, VT>::RBasis(MT& A, MT& u) {
	cout << "Start computing R basis..." << endl;
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
void SPH_RIDG<pT, MT, VT>::QBasis(MT& Q, MT& u) {
	cout << "Start computing Q basis..." << endl;
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
void SPH_RIDG<pT, MT, VT>::normBasis(MT& mat) {
	MT e = mat * mat.transpose();
	SelfAdjointEigenSolver<MT> eigensolver(mat.rows());
	eigensolver.compute(e, EigenvaluesOnly);
	pT lVal = eigensolver.eigenvalues()[mat.rows() - 1];
	mat = (1 / sqrt(lVal)) * mat;
}

#endif