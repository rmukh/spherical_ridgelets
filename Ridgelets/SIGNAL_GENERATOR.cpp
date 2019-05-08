#include "SIGNAL_GENERATOR.h"

SIGNAL_GENERATOR::SIGNAL_GENERATOR() : inputVolume(NULL), nGradImgs(0), nOfImgs(0) {}

SIGNAL_GENERATOR::SIGNAL_GENERATOR(string & iv) : inputVolume(iv), nGradImgs(0), nOfImgs(0) {}

SIGNAL_GENERATOR::~SIGNAL_GENERATOR() {}

bool SIGNAL_GENERATOR::is_path_exists(const string &s)
{
	struct stat buffer;
	return (stat(s.c_str(), &buffer) == 0);
}

int SIGNAL_GENERATOR::readVolume(MatrixType & GradientDirections, DiffusionImagePointer & image) {
	GradientDirections.resize(0, 3);
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
			cout << "NDHR or NRRD file formats only! Please, check file type." << endl;
			return EXIT_FAILURE;
		}
	}
	if (!is_path_exists(inputVolume)) {
		cout << "Input dMRI image is not exists! Please, check the path and file name." << endl;
		return EXIT_FAILURE;
	}

	// Get image
	reader->SetFileName(inputVolume);
	try
	{
		reader->Update();
		image = reader->GetOutput();
	}
	catch (itk::ExceptionObject)
	{
		cerr << "Can't read input dMRI file! Please, check that file is not corrupted." << endl;
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
			if (x == 0.0 && y == 0.0 && z == 0.0)
				continue;

			GradientDirections.conservativeResize(GradientDirections.rows() + 1, GradientDirections.cols());
			GradientDirections.row(GradientDirections.rows() - 1) << x, y, z;
			++nGradImgs;
		}
		else if (itKey->find("DWMRI_b-value") != string::npos)
		{
			//cout << *itKey << " -> " << metaString << endl;
			is_b0 = true;
			b0 = stod(metaString.c_str());
		}
	}
	// Normalize gradients
	GradientDirections.rowwise().normalize();

	itk::Size<3> img_size = image->GetLargestPossibleRegion().GetSize();
	cout << "Image size: (" << img_size[0] << " " << img_size[1] << " " << img_size[2] << "). ";
	cout << "Number of gradient images: " << nGradImgs << ". Number of reference images: " << nOfImgs - nGradImgs << endl;
	cout << "b-value " << b0 << endl;

	if (!is_b0)
	{
		cerr << "b-value not specified in file's header." << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int SIGNAL_GENERATOR::ExtractMatrix(MaskImagePointer &mask, MatrixType &signal, MatrixType &grad_dirs)
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
	precisionType avrg_b0 = -1.0;
	unsigned first_grad_image_index = nOfImgs - nGradImgs;
	if (first_grad_image_index == 0) {
		cout << "Warning! There is no baseline image in the input file!" << endl;
		avrg_b0 = 0.0;
	}
	signal = MatrixType::Zero(nGradImgs, N_of_voxels);
	DiffusionImageType::PixelType voxel_content;

	// Iterate over all voxels
	if (mask) {
		ConstIterator it_i(img, img->GetRequestedRegion());
		MaskIterator it_m(mask, mask->GetRequestedRegion());
		unsigned vox = 0;

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

					// Calculate average b0 value for the voxel
					if (avrg_b0 == 0.0)
					{
						avrg_b0 = 1.0;
					}
					else
					{
						precisionType avrgb0_sum = 0;
						for (unsigned i = 0; i < first_grad_image_index; ++i) {
							precisionType vol = voxel_content.GetElement(i);
							avrgb0_sum += vol;
						}
						avrg_b0 = avrgb0_sum / first_grad_image_index;
					}

					for (unsigned i = first_grad_image_index; i < nOfImgs; ++i) {
						// Get diffusion pixel and normalize to b0
						precisionType vol;
						if (avrg_b0 == 0.0)
							vol = 1e-10;
						else
							vol = voxel_content.GetElement(i) / avrg_b0;
						// Remove negative values
						if (vol > 0)
							signal(i - first_grad_image_index, vox) = vol;
					}
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
		unsigned vox = 0;
		it.SetDirection(0);
		it.GoToBegin();
		while (!it.IsAtEnd())
		{
			while (!it.IsAtEndOfLine())
			{
				voxel_content = it.Get();

				// Calculate average b0 value for the voxel
				if (avrg_b0 == 0.0)
				{
					avrg_b0 = 1.0;
				}
				else
				{
					precisionType avrgb0_sum = 0;
					for (unsigned i = 0; i < first_grad_image_index; ++i) {
						precisionType vol = voxel_content.GetElement(i);
						avrgb0_sum += vol;
					}
					avrg_b0 = avrgb0_sum / first_grad_image_index;
				}

				for (unsigned i = first_grad_image_index; i < nOfImgs; ++i) {
					// Get diffusion pixel and normalize to b0
					precisionType vol;
					if (avrg_b0 == 0.0)
						vol = 1e-10;
					else
						vol = voxel_content.GetElement(i) / avrg_b0;
					// Remove negative values
					if (vol > 0)
						signal(i - first_grad_image_index, vox) = vol;
				}

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
