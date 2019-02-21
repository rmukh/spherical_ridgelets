#ifndef SPH_RIDG_H
#define SPH_RIDG_H

#include "UtilMath.h"

class SPH_RIDG
{
public:
	/* 
	For more information on spherical ridgelets, please refere to:
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3073602/
	*/
	int J; //J value
	double rho; //rho value
	UtilMath UM = UtilMath(); //Object of UtilMath class
	MatrixType C;

	SPH_RIDG();
	SPH_RIDG(unsigned JIn, double rhoIn);
	void init(); //pre-compute all necessary matricies/vectors

	MatrixType RBasis(MatrixType u); //return spherical ridgelets basis matrix
	MatrixType QBasis(MatrixType u); //For visualizing pupose only
	MatrixType normBasis(MatrixType B); //normalize basis

private:
	int mcut;
	MatrixType h;
	MatrixType psi;
	MatrixType t;
	MatrixType Lmd;
	double tau;
	int m0;
	VectorXi M0;
};

#endif // !SPH_RIDG_H



