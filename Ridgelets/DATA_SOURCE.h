#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include "UtilMath.h"
#include "SPH_RIDG.h"

#include <string>
#include <fstream>

class DATA_SOURCE
{
public:
	int readNRRD(string inputVolume, MatrixType &GradientDirections, DiffusionImagePointer &image, unsigned &nGradImgs, unsigned &nOfImgs);
	void testFNC();
	void readTestData(MatrixType& g, MatrixType& s);
	void DWI2Matrix(DiffusionImagePointer & img, MatrixType & signal, unsigned & nGradImgs, unsigned & nOfImgs);
	void matrixToFile(const string & fname, MatrixType & matrix);
	void fileToMatrix(const string & fname, MatrixType & arr);
	template<typename T>
	void printVec(const string & name, vector<T>& v);
};

#endif
