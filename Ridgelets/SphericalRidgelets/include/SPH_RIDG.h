#ifndef SPH_RIDG_H
#define SPH_RIDG_H

#include "UtilMath.h"

//precisionType - pT, MatrixType - MT
template <class pT, class MT, class VT>
class SPH_RIDG
{
public:
	/* 
	For more information on spherical ridgelets, please refere to:
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073602/
	*/
	int J; //J value
	pT rho; //rho value
	UtilMath<pT, MT, VT> UM; //Object of UtilMath class
	MT C;

	SPH_RIDG();
	~SPH_RIDG();

	SPH_RIDG(unsigned JIn, pT rhoIn);
	void init(); //pre-compute all necessary matricies/vectors

	void RBasis(MT& A, MT& u); //return spherical ridgelets basis matrix
	void QBasis(MT& Q, MT& u); //For visualizing pupose only
	const MT& getScaleWeights() const;
	void FindMaxRidgeletMaxInDMRI(MT& fin, MT& R, pT nms_angle_degrees = static_cast<pT>(20.0));
	bool TestDirectRidgeletMaxima(pT& abs_dot, pT& score, int& coeff_index);
	void normBasis(MT& mat); //normalize basis
	const pT USER_RHO_DEFAULT = static_cast<pT>(3.125);
private:
	const pT TAIL_THRESHOLD = static_cast<pT>(1e-6);
	int mcut;
	MT h;
	MT psi;
	MT scale_weights;
	MT t;
	MT Lmd;
	pT tau;
	int m0;
	VectorXi M0;
};

#include "SPH_RIDG.hpp"

#endif // !SPH_RIDG_H

