#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include "UtilMath.h"
#include "SPH_RIDG.h"
#include "SIGNAL_GENERATOR.h"

class DATA_SOURCE
{
public:
	struct input_parse {
		string input_dmri, input_mask, output_ridgelets, output_fiber_max_odf, output_odf;
		unsigned lvl;
		bool is_compress;
	};

	DATA_SOURCE();
	~DATA_SOURCE();

	int CLI(int argc, char * argv[], input_parse & output);
	int readMask(string inputMask, MaskImagePointer & image);
	void set_header(DiffusionImagePointer & dest);
	bool is_path_exists(const string & s);
	void testFNC();
	void readTestData(MatrixType & g, MatrixType & s);
	int DWI2Matrix(string & dmri_file, MaskImagePointer & mask, MatrixType & signal, MatrixType & grad_dirs);
	void Matrix2DWI(DiffusionImagePointer & img, MaskImagePointer & mask, MatrixType & arr);
	void matrixToFile(const string & fname, MatrixType & matrix);
	void fileToMatrix(const string & fname, MatrixType & arr);
	template<typename D>
	int save_to_file(const string & fname, typename D::Pointer & image, bool is_compress);
	template<typename T>
	void printVec(const string & name, vector<T>& v);
	dMRI_h_info h;
};

#endif
