#ifndef UTILMATH_H
#define UTILMATH_H

#include <iostream>
#include "rdgls_types.h"

using namespace std;

class UtilMath {
public:
	//Useful constants
	const double PI = atan(1.0) * 4.0;

	//Useful pre-defined types

	//Functions
	void spiralsample(MatrixType& u, unsigned flg, unsigned N);
	void fura(MatrixType& Lmd, unsigned n);
	void polyleg(MatrixType& P, MatrixType x, unsigned n);
	MatrixType convhulln(MatrixType & u);
};

#endif // ! UTILMATH_H