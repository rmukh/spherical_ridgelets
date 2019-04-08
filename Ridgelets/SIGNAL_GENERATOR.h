#ifndef SIGNAL_GENERATOR_H
#define SIGNAL_GENERATOR_H

#include "rdgls_types.h"

class SIGNAL_GENERATOR
{
public:
	SIGNAL_GENERATOR();
	SIGNAL_GENERATOR(string & iv);
	~SIGNAL_GENERATOR();

	bool is_path_exists(const string & s);
	int readVolume(MatrixType & GradientDirections, DiffusionImagePointer & image);
	int ExtractMatrix(MaskImagePointer & mask, MatrixType & signal, MatrixType & grad_dirs);
	unsigned nGradImgs = 0; // Number of gradient images
	unsigned nOfImgs = 0; // Total number of images (including b0)
	dMRI_h_info h;

private:
	const string & inputVolume;
};

#endif // !SIGNAL_GENERATOR

