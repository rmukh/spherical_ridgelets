#include "SIGNAL_GENERATOR.h"

SIGNAL_GENERATOR::SIGNAL_GENERATOR(string& iv) : nGradImgs(0), nOfImgs(0), inputVolume(iv) {}

bool SIGNAL_GENERATOR::is_path_exists(const string& s)
{
	struct stat buffer;
	return (stat(s.c_str(), &buffer) == 0);
}

precisionType SIGNAL_GENERATOR::calculate_average_b0(DiffusionImageType::PixelType voxel_content, unsigned first_grad_image_index) {
	// Calculate average b0 value for the voxel

	precisionType average_b0_sum = sphZero;
	precisionType average_b0 = sphOne;
	precisionType volume_element = sphZero;

	for (unsigned i = 0; i < first_grad_image_index; ++i) {
		volume_element = voxel_content.GetElement(i);
		average_b0_sum += volume_element;
	}
	average_b0 = average_b0_sum / first_grad_image_index;

	return average_b0;
}

void SIGNAL_GENERATOR::compute_n_fill_voxel_element(
	MatrixType& signal,
	DiffusionImageType::PixelType voxel_content,
	unsigned first_grad_image_index,
	precisionType average_b0,
	unsigned vox,
	bool no_b0
) {
	// iterate over diffusion-encoding pixels

	precisionType volume_element = sphZero;
	precisionType vol_elem_eps = 1e-10;

	for (unsigned i = first_grad_image_index; i < nOfImgs; ++i) {
		// Get diffusion pixel and normalize to b0
		volume_element = voxel_content.GetElement(i);
		if (fabs(average_b0 - sphZero) <= vol_elem_eps)
			volume_element = vol_elem_eps;
		else if (!no_b0)
			volume_element /= average_b0;
		// Remove negative values
		if (volume_element > 0)
			signal(i - first_grad_image_index, vox) = volume_element;
	}
}

int SIGNAL_GENERATOR::readVolume(MatrixType& GradientDirections, DiffusionImagePointer& image) {
	bool ext_grad_use = true;
	if (GradientDirections.size() == 0) {
		GradientDirections.resize(0, 3);
		ext_grad_use = false;
	}

	bool is_b0 = false;
	precisionType b0 = 0;
	precisionType x, y, z;

	// Temporary register factories cause don't use cmake
	itk::NrrdImageIOFactory::RegisterOneFactory();

	// Another way to store image data
	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->ReleaseDataFlagOn();

	// Make some inputfiles checks
	string ext_vol = inputVolume.substr(inputVolume.length() - 4, inputVolume.length());

	if (ext_vol.compare("nhdr")) {
		if (ext_vol.compare("nrrd")) {
			cerr << "NDHR or NRRD file formats only! Please, check file type." << endl;
			return EXIT_FAILURE;
		}
	}
	if (!is_path_exists(inputVolume)) {
		cerr << "Input dMRI image is not exists! Please, check the path and file name." << endl;
		return EXIT_FAILURE;
	}

	// Get image
	reader->SetFileName(inputVolume);
	try
	{
		reader->Update();
		image = reader->GetOutput();
	}
	catch (itk::ExceptionObject& ex)
	{
		cerr << "Can't read input dMRI file! Please, check that file is not corrupted." << endl;
		cerr << "Extra error message: " << ex << endl;
		return EXIT_FAILURE;
	}

	// Get and process header information
	h.spc_h = image->GetSpacing();
	h.dirs_h = image->GetDirection();
	h.orig_h = image->GetOrigin();
	h.reg_h = image->GetLargestPossibleRegion();
	h.comp_h = image->GetNumberOfComponentsPerPixel();

	itk::MetaDataDictionary imgMetaDictionary = image->GetMetaDataDictionary();
	vector<string> imgMetaKeys = imgMetaDictionary.GetKeys();
	vector<string>::const_iterator itKey = imgMetaKeys.begin();
	string metaString;

	// Get diffusion-encoding direction from dmri volume if an external file is not provided
	if (!ext_grad_use) {
		for (; itKey != imgMetaKeys.end(); ++itKey)
		{
			itk::ExposeMetaData<string>(imgMetaDictionary, *itKey, metaString);
			if (itKey->find("DWMRI_gradient") != string::npos)
			{
				//cout << *itKey << " -> " << metaString << endl;
				istringstream ss(metaString.c_str());
				ss >> x; ss >> y; ss >> z;

				++nOfImgs;
				// If the direction is 0.0, this is a reference image
				if (x == sphZero && y == sphZero && z == sphZero)
					continue;

				GradientDirections.conservativeResize(GradientDirections.rows() + 1, GradientDirections.cols());
				GradientDirections.row(GradientDirections.rows() - 1) << x, y, z;
				++nGradImgs;
			}
		}
	}
	else {
		for (int i = 0; i < GradientDirections.rows(); ++i) {
			++nOfImgs;
			if (GradientDirections(i, 0) == sphZero && GradientDirections(i, 1) == sphZero && GradientDirections(i, 2) == sphZero)
				continue;
			++nGradImgs;
		}
		for (int i = 0; i < GradientDirections.rows(); ++i) {
			if (GradientDirections(i, 0) == sphZero && GradientDirections(i, 1) == sphZero && GradientDirections(i, 2) == sphZero)
				remove_row(GradientDirections, i);
		}
	}

	// Get b-value from the dmri volume
	itKey = imgMetaKeys.begin(); // reset iterator's pointer
	for (; itKey != imgMetaKeys.end(); ++itKey)
	{
		itk::ExposeMetaData<string>(imgMetaDictionary, *itKey, metaString);
		if (itKey->find("DWMRI_b-value") != string::npos)
		{
			// cout << *itKey << " -> " << metaString << endl;
			is_b0 = true;
			b0 = stod(metaString.c_str());
		}
	}

	// Normalize gradients
	GradientDirections.rowwise().normalize();

	itk::Size<3> img_size = image->GetLargestPossibleRegion().GetSize();
	unsigned short int n_refs = nOfImgs - nGradImgs;
	cout << "Image size: (" << img_size[0] << " " << img_size[1] << " " << img_size[2] << "). ";
	cout << "Number of gradient images: " << nGradImgs << ". Number of reference images: " << n_refs << endl;
	if (!is_b0) {
		cout << "Be cautious! b-value is not specified in file's header. Sometimes it indicates that the input is corrupted." << endl;
	}
	else {
		cout << "b-value " << b0 << endl;
	}

	return EXIT_SUCCESS;
}

int SIGNAL_GENERATOR::ExtractMatrix(MaskImagePointer& mask, MatrixType& signal, MatrixType& grad_dirs)
{
	DiffusionImagePointer img = DiffusionImageType::New();
	int res_dmri = readVolume(grad_dirs, img);
	if (res_dmri)
		return EXIT_FAILURE;

	int N_of_voxels = 0;

	if (mask) {
		//If mask provided

		MaskIterator m_it(mask, mask->GetRequestedRegion());

		m_it.SetDirection(0);
		m_it.GoToBegin();
		while (!m_it.IsAtEnd())
		{
			while (!m_it.IsAtEndOfLine())
			{
				if (m_it.Get() == 1) {
					N_of_voxels += 1;
				}
				++m_it;
			}
			m_it.NextLine();
		}
	}
	else {
		DiffusionImageType::SizeType sz = img->GetLargestPossibleRegion().GetSize();
		N_of_voxels = sz[0] * sz[1] * sz[2];
	}

	// some in-dwi-loop variables declaration
	precisionType vol_elem = sphZero;
	precisionType avrg_b0 = sphOne;
	bool no_b0 = false;
	unsigned vox = 0;

	unsigned first_grad_image_index = nOfImgs - nGradImgs;
	if (first_grad_image_index == 0) {
		cout << "Warning! There is no baseline image in the input file!" << endl;
		cout << "Make sure that all diffusion-encoded volumes were normalized by an average b0!" << endl;
		no_b0 = true;
	}
	signal = MatrixType::Zero(nGradImgs, N_of_voxels);
	DiffusionImageType::PixelType voxel_content;

	// Iterate over all voxels
	if (mask) {
		ConstIterator it_i(img, img->GetRequestedRegion());
		MaskIterator it_m(mask, mask->GetRequestedRegion());

		it_i.SetDirection(0);
		it_i.GoToBegin();
		it_m.SetDirection(0);
		it_m.GoToBegin();

		while (!it_i.IsAtEnd())
		{
			while (!it_i.IsAtEndOfLine())
			{
				// Take pixels only within the mask region
				if (it_m.Get() == 1) {
					voxel_content = it_i.Get();

					if (!no_b0)
						avrg_b0 = calculate_average_b0(voxel_content, first_grad_image_index);

					compute_n_fill_voxel_element(signal, voxel_content, first_grad_image_index, avrg_b0, vox, no_b0);

					++vox;
				}
				++it_i;
				++it_m;
			}
			it_i.NextLine();
			it_m.NextLine();
		}
	}
	else {
		ConstIterator it(img, img->GetRequestedRegion());
		it.SetDirection(0);
		it.GoToBegin();
		while (!it.IsAtEnd())
		{
			while (!it.IsAtEndOfLine())
			{
				voxel_content = it.Get();

				if (!no_b0)
					avrg_b0 = calculate_average_b0(voxel_content, first_grad_image_index);

				compute_n_fill_voxel_element(signal, voxel_content, first_grad_image_index, avrg_b0, vox, no_b0);

				++vox;
				++it;
			}
			it.NextLine();
		}
	}

	cout << "Total number of voxels to process (with mask): " << signal.cols() << endl << endl;
	img = NULL;
	return EXIT_SUCCESS;
}

void SIGNAL_GENERATOR::remove_row(MatrixType& a, Eigen::Index del)
{
	unsigned cols = a.cols();
	unsigned rows = a.rows() - 1;

	if (del < rows)
		a.block(del, 0, rows - del, cols) = a.block(del + 1, 0, rows - del, cols);

	a.conservativeResize(rows, cols);
}
