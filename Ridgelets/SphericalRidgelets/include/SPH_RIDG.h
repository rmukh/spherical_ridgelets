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
	UtilMath<pT, MT> UM; //Object of UtilMath class
	MT C;

	SPH_RIDG();
	~SPH_RIDG();

	SPH_RIDG(unsigned JIn, pT rhoIn);
	void init(); //pre-compute all necessary matricies/vectors

	void RBasis(MT& A, MT& u); //return spherical ridgelets basis matrix
	void QBasis(MT& Q, MT& u); //For visualizing pupose only
	void normBasis(MT& mat); //normalize basis

private:
	int mcut;
	MT h;
	MT psi;
	MT t;
	MT Lmd;
	pT tau;
	int m0;
	VectorXi M0;
};

#include "../src/SPH_RIDG.cpp"

#endif // !SPH_RIDG_H



