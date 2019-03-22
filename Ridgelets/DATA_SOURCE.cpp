#include "DATA_SOURCE.h"

DATA_SOURCE::DATA_SOURCE() {}
DATA_SOURCE::~DATA_SOURCE() { cout << "DATA_SOURCE destructed" << endl; }

int DATA_SOURCE::CLI(int argc, char* argv[], input_parse& output) {
	if (argc < 5)
	{
		cerr << "Usage: Ridgelets -i dMRI_file and at least one output: -ridg, -odf, -omd" << endl;
		cerr << "Optional input arguments: -m mask_file" << endl;
		cerr << "Possible output argumet(s): -ridg ridgelet_file -odf ODF_values -omd ODF_maxima_dir_&_value -c enable compression" << endl;
		return EXIT_FAILURE;
	}

	bool inp1 = false;
	bool out1 = false;
	output.is_compress = false;
	output.lvl = 4;
	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "-i")) {
			output.input_dmri = argv[i + 1];
			inp1 = true;
		}
		if (!strcmp(argv[i], "-m")) {
			output.input_mask = argv[i + 1];
		}
		if (!strcmp(argv[i], "-lvl")) {
			float order = stof(argv[i + 1]);
			if (order == floor(order) && order > 0) {
				output.lvl = order;
			}
			else {
				cout << "The value for icosahedron "
					"tesselation order provided is in the wrong "
					"format (must be a positive integer). "
					"So default value 4 used." << endl;
			}
		}
		if (!strcmp(argv[i], "-ridg")) {
			output.output_ridgelets = argv[i + 1];
			out1 = true;
		}
		if (!strcmp(argv[i], "-odf")) {
			output.output_odf = argv[i + 1];
			out1 = true;
		}
		if (!strcmp(argv[i], "-omd")) {
			output.output_fiber_max_odf = argv[i + 1];
			out1 = true;
		}
		if (!strcmp(argv[i], "-c")) {
			output.is_compress = true;
		}
	}
	if (!inp1 || !out1) {
		cerr << "Please, provide at least one input AND one output file names" << endl;
		return EXIT_FAILURE;
	}
	return 0;
}

int DATA_SOURCE::readMask(string inputMask, MaskImagePointer& image) {
	// We need mask within the main program so it is implemented in that class

	// Temporary register factories cause don't use cmake
	itk::NrrdImageIOFactory::RegisterOneFactory();

	// Another way to store image data
	MaskReaderType::Pointer reader = MaskReaderType::New();

	// Make some inputfiles checks
	if (!is_path_exists(inputMask)) {
		cout << "Input mask image is not provided. Please, stop program and provide mask file if you forget to include it." << endl;
		return EXIT_SUCCESS;
	}

	string ext_vol = inputMask.substr(inputMask.length() - 4, inputMask.length());

	if (ext_vol.compare("nhdr")) {
		if (ext_vol.compare("nrrd")) {
			cout << "NDHR or NRRD file formats only! Please, check file type." << endl;
			return EXIT_FAILURE;
		}
	}

	// Get image
	reader->SetFileName(inputMask);
	try
	{
		reader->Update();
		image = reader->GetOutput();
	}
	catch (itk::ExceptionObject)
	{
		cerr << "Can't read input mask file! Please, check that file is not corrupted." << endl;
		return EXIT_FAILURE;
	}
	cout << "Mask file loaded" << endl;
	return EXIT_SUCCESS;
}

template<typename D>
int DATA_SOURCE::save_to_file(const string& fname, typename D::Pointer& image, bool is_compress) {
	typedef itk::ImageFileWriter<D> ImageWriterType;
	typename ImageWriterType::Pointer writer = ImageWriterType::New();

	writer->SetFileName(fname);
	writer->SetInput(image);
	writer->SetUseCompression(is_compress);

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		cerr << "Error while saving file on disk:\n " << err << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

void DATA_SOURCE::set_header(DiffusionImagePointer& dest) {
	//dest->SetMetaDataDictionary(src->GetMetaDataDictionary());
	//dest->SetSpacing(src->GetSpacing());
	//dest->SetDirection(src->GetDirection());
	//dest->SetOrigin(src->GetOrigin());
	//dest->SetRegions(src->GetLargestPossibleRegion());
	//dest->SetNumberOfComponentsPerPixel(src->GetNumberOfComponentsPerPixel());
}

bool DATA_SOURCE::is_path_exists(const string &s)
{
	struct stat buffer;
	return (stat(s.c_str(), &buffer) == 0);
}

void DATA_SOURCE::testFNC() {
	//MatrixType test(19, 1);
	//test << -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;
	//MatrixType res = UM.polyleg(test, 5);
	SPH_RIDG test;
	//cout << test.M0;
	UtilMath UM = UtilMath();
	MatrixType g;
	UM.spiralsample(g, 2, 100);
	MatrixType A;
	test.RBasis(A, g);
	test.normBasis(A);

	test.normBasis(A);
}

void DATA_SOURCE::readTestData(MatrixType& g, MatrixType& s) {
	// Works only on my local PC!
	fileToMatrix("D:\\test\\bvec.txt", g);
	fileToMatrix("D:\\test\\sign.txt", s);

	//normalize
	MatrixType gnorm = g.array().pow(2).rowwise().sum().sqrt();
	for (int i = 0; i < 51; ++i)
		g.row(i) = g.row(i) / gnorm(i);
}

int DATA_SOURCE::DWI2Matrix(string &dmri_file, MaskImagePointer &mask, MatrixType &signal, MatrixType &grad_dirs)
{
	SIGNAL_GENERATOR sg(dmri_file);
	int res_dmri = sg.ExtractMatrix(mask, signal, grad_dirs);
	if (res_dmri)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}

void DATA_SOURCE::Matrix2DWI(DiffusionImagePointer &img, MaskImagePointer &mask, MatrixType &arr) {
	unsigned n_of_components = arr.rows();

	// Make a vector of zeros
	VariableVectorType zeros_vec;
	zeros_vec.SetSize(n_of_components);
	for (unsigned i = 0; i < n_of_components; ++i) {
		zeros_vec[i] = 0;
	}
	img->FillBuffer(zeros_vec);

	VariableVectorType vec_to_fill;
	vec_to_fill.SetSize(n_of_components);

	Iterator it(img, img->GetRequestedRegion());

	it.SetDirection(0);
	it.GoToBegin();
	unsigned vox = 0;

	if (mask != nullptr) {
		// Mask iterator
		MaskIterator it_m(mask, mask->GetRequestedRegion());

		it_m.SetDirection(0);
		it_m.GoToBegin();

		while (!it.IsAtEnd())
		{
			while (!it.IsAtEndOfLine())
			{
				if (it_m.Get() == 1) {
					for (unsigned i = 0; i < n_of_components; ++i) {
						vec_to_fill[i] = arr(i, vox);
					}
					it.Set(vec_to_fill);
					++vox;
				}
				++it;
				++it_m;
			}
			it.NextLine();
			it_m.NextLine();
		}
	}
	else {
		while (!it.IsAtEnd())
		{
			while (!it.IsAtEndOfLine())
			{
				for (unsigned i = 0; i < n_of_components; ++i) {
					vec_to_fill[i] = arr(i, vox);
				}

				it.Set(vec_to_fill);
				++it;
				++vox;
			}
			it.NextLine();
		}
	}
}

void DATA_SOURCE::matrixToFile(const string& fname, MatrixType& matrix) {
	/*
	Write space separated matrix of MatrixType to text type file. Useful to debugging.
	*/
	ofstream file(fname);
	if (file.is_open())
		file << matrix << '\n';
}

void DATA_SOURCE::fileToMatrix(const string& fname, MatrixType& matrix)
{
	/*
	Read space separated matricies and vectors. Useful to debug with matlab code.
	*/
	unsigned number_of_rows = 0;
	unsigned number_of_cols = 0;
	std::string line;

	ifstream infile;
	infile.open(fname);

	// Get number of lines
	while (getline(infile, line))
		++number_of_rows;

	// Get number of columns
	infile.clear();
	infile.seekg(0, ios::beg);
	getline(infile, line);
	stringstream stream(line);
	while (stream)
	{
		std::string c;
		stream >> c;
		if (c.length())
			++number_of_cols;
	}

	matrix.resize(number_of_rows, number_of_cols);

	// Filling matrix
	infile.clear();
	infile.seekg(0, ios::beg);

	for (unsigned i = 0; i < number_of_rows; ++i)
	{
		getline(infile, line);
		stringstream stream(line);

		for (unsigned j = 0; j < number_of_cols; ++j)
			stream >> matrix(i, j);
	}
}

template <typename T>
void DATA_SOURCE::printVec(const string& name, vector<T>& v) {
	cout << name << endl;
	for (auto i : v)
		std::cout << i << ' ';
	cout << endl;
}

// Necessary templates for utilized types. Needed because implementation and declaration of function in separate files
template void DATA_SOURCE::printVec<int>(const string&, vector<int>&);
template void DATA_SOURCE::printVec<unsigned>(const string&, vector<unsigned>&);
template void DATA_SOURCE::printVec<Eigen::Index>(const string&, vector<Eigen::Index>&);

template int DATA_SOURCE::save_to_file<DiffusionImageType>(const string&, DiffusionImageType::Pointer&, bool);
