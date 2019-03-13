#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include "UtilMath.h"
#include "SPH_RIDG.h"

#include <string>
#include <fstream>

class DATA_SOURCE
{
public:
	int readVolume(string inputVolume, MatrixType & GradientDirections, DiffusionImagePointer & image, unsigned & nGradImgs, unsigned & nOfImgs);
	int readMask(string inputMask, MaskImagePointer & image);
	void copy_header(DiffusionImagePointer & src, DiffusionImagePointer & dest);
	bool is_path_exists(const string & s);
	void testFNC();
	void readTestData(MatrixType & g, MatrixType & s);
	void DWI2Matrix(DiffusionImagePointer & img, MaskImagePointer & mask, MatrixType & signal, unsigned & nGradImgs, unsigned & nOfImgs);
	void Matrix2DWI(DiffusionImagePointer & img, MaskImagePointer & mask, MatrixType & arr);
	void matrixToFile(const string & fname, MatrixType & matrix);
	void fileToMatrix(const string & fname, MatrixType & arr);
	template<typename D>
	int save_to_file(const string & fname, typename D::Pointer & image, bool is_compress);
	template<typename T>
	void printVec(const string & name, vector<T>& v);
};

#endif
