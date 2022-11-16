#ifndef SIGNAL_GENERATOR_H
#define SIGNAL_GENERATOR_H

#include "rdgls_types.h"

class SIGNAL_GENERATOR
{
public:
	SIGNAL_GENERATOR(string& iv);

	bool is_path_exists(const string& s);
	precisionType calculate_average_b0(DiffusionImageType::PixelType voxel_content, unsigned first_grad_image_index);
	void compute_n_fill_voxel_element(
		MatrixType& signal,
		DiffusionImageType::PixelType voxel_content,
		unsigned first_grad_image_index,
		precisionType average_b0,
		unsigned vox,
		bool no_b0
	);
	int readVolume(MatrixType& GradientDirections, DiffusionImagePointer& image);
	int ExtractMatrix(MaskImagePointer& mask, MatrixType& signal, MatrixType& grad_dirs);
	void remove_row(MatrixType& a, Eigen::Index del);
	unsigned nGradImgs; // Number of gradient images
	unsigned nOfImgs; // Total number of images (including b0)
	dMRI_h_info h;

private:
	const string& inputVolume;
};

#endif // !SIGNAL_GENERATOR

