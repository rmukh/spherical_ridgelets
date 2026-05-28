#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include "rdgls_types.h"
#include "SPH_RIDG.h"
#include "SIGNAL_GENERATOR.h"

class DATA_SOURCE
{
public:
	typedef struct input_parse {
		string input_dmri, input_mask, output_ridgelets, output_fiber_max_odf,
			output_fiber_max_ridgelets, output_odf, signal_recon, output_A,
			external_gradients, ext_signal_recon;
		precisionType max_odf_thresh = static_cast<precisionType>(0.7);
		precisionType fista_lambda = static_cast<precisionType>(0.01);
		precisionType sph_rho = static_cast<precisionType>(3.125);
		precisionType ridgelet_nms_angle = static_cast<precisionType>(20.0);
		unsigned lvl = 4;
		unsigned sph_J = 2;
		int n_splits = -1;
		int nth = -1;
		bool is_compress = false;
		bool print_scale_weights = false;
		bool test_direct_ridgelet_maxima = false;
		int fista_iterations = 2000;
		precisionType fista_tolerance = static_cast<precisionType>(0.001);
	} input_parse;

	DATA_SOURCE();
	~DATA_SOURCE();

	int CLI(int argc, char * argv[], input_parse * output);
	void short_summary(const input_parse & params);
	int readMask(string inputMask, MaskImagePointer & image);
	void set_header(DiffusionImagePointer & dest);
	bool is_path_exists(const string & s);
	void testFNC();
	void readTestData(MatrixType & g, MatrixType & s);
	void data_saving_info_out(unsigned long int coef_size, string name);
	void estimate_memory(MatrixType & s, MatrixType & A, const input_parse & params);
	int compute_splits(unsigned s_size);
	int DWI2Matrix(input_parse* input_args, MaskImagePointer & mask, MatrixType & signal, MatrixType & grad_dirs);
	void Matrix2DWI(DiffusionImagePointer & img, MaskImagePointer & mask, MatrixType & arr);
	void matrixToFile(const string & fname, MatrixType & matrix);
	void fileToMatrix(const string & fname, MatrixType & arr);
	void fileGradientsToMatrix(const string& fname, MatrixType& arr);
	template<typename D>
	int save_to_file(const string & fname, typename D::Pointer & image, bool is_compress);
	template<typename T>
	void printVec(const string & name, vector<T>& v);
	dMRI_h_info h;
};

#endif
