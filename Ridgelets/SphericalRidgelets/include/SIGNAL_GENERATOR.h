#ifndef SIGNAL_GENERATOR_H
#define SIGNAL_GENERATOR_H

#include "rdgls_types.h"

class SIGNAL_GENERATOR
{
public:
	SIGNAL_GENERATOR(string& iv);

	bool is_path_exists(const string& s);
	string encoding_direction_2_string(const MatrixType& a);
	precisionType calculate_average_b0(DiffusionImageType::PixelType voxel_content);
	void compute_n_fill_voxel_element(
		MatrixType& signal,
		DiffusionImageType::PixelType voxel_content,
		precisionType average_b0,
		unsigned vox,
		bool no_b0
	);
	int readVolume(MatrixType& GradientDirections, DiffusionImagePointer& image);
	int ExtractMatrix(MaskImagePointer& mask, MatrixType& signal, MatrixType& grad_dirs);
	void remove_row(MatrixType& a, Eigen::Index del);
	vector<unsigned> gradImgs; // gradient images indices
	vector<unsigned> b0Imgs; // b0 images indices
	dMRI_h_info h;

private:
	const string& inputVolume;
};

#endif // !SIGNAL_GENERATOR

