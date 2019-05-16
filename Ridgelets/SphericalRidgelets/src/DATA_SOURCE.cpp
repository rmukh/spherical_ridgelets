#include "DATA_SOURCE.h"

DATA_SOURCE::DATA_SOURCE() {}
DATA_SOURCE::~DATA_SOURCE() {}

int DATA_SOURCE::CLI(int argc, char* argv[], input_parse& output) {
	if (argc < 5)
	{
		cerr << "Usage: Ridgelets -i dMRI file AND at least one output: -ridg, -sr, -odf, -omd" << endl;
		cerr << "Optional input arguments: -m mask file, -lvl ridgelets order, -nspl splits "
			"coefficient, -mth maxima ODF threshold, -lmd FISTA lambda, -sj Spherical ridgelets J, "
			"-srho Spherical ridgelets rho, -nth number of threads to use" << endl;
		cerr << "Possible output argumet(s): -ridg ridgelet_file, -sr signal reconstruction, -odf ODF_values, -omd ODF_maxima_dir_&_value, -c enable compression" << endl;
		return EXIT_FAILURE;
	}

	bool inp1 = false;
	bool out1 = false;
	output.is_compress = false;
	output.lvl = 4;
	output.n_splits = -1;
	output.max_odf_thresh = 0.7;
	output.fista_lambda = 0.01;
	output.sph_J = 2;
	output.sph_rho = 3.125;
	output.nth = -1;
	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "-i")) {
			output.input_dmri = argv[i + 1];
			inp1 = true;
		}
		if (!strcmp(argv[i], "-m")) {
			output.input_mask = argv[i + 1];
		}
		if (!strcmp(argv[i], "-sj")) {
			float sj = stof(argv[i + 1]);
			if (sj == floor(sj) && sj > 0) {
				output.sph_J = sj;
			}
			else {
				cout << "The J value of spherical ridgelets  "
					"basis provided is in the wrong "
					"format (must be a positive integer). "
					"So, default value 2 used." << endl;
			}
		}
		if (!strcmp(argv[i], "-srho")) {
			float rho = stof(argv[i + 1]);
			if (rho > 0) {
				output.sph_rho = rho;
			}
			else {
				cout << "The maxima ODF search threshold "
					"coefficient provided is in the wrong "
					"format (must be float point number > 0). "
					"So, the default value 3.125 used." << endl;
			}
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
					"So, default value 4 used." << endl;
			}
		}
		if (!strcmp(argv[i], "-nspl")) {
			float splt = stof(argv[i + 1]);
			if (splt == floor(splt) && splt > 0) {
				output.n_splits = splt;
			}
			else {
				cout << "The split coefficient for the ridgelets "
					"computation provided is in the wrong "
					"format (must be a positive integer). "
					"So, the default value (available CPU threads * 2) used." << endl;
			}
		}
		if (!strcmp(argv[i], "-mth")) {
			float th = stof(argv[i + 1]);
			if (th > 0 && th < 1) {
				output.max_odf_thresh = th;
			}
			else {
				cout << "The maxima ODF search threshold "
					"coefficient provided is in the wrong "
					"format (must be in (0, 1)). "
					"So, the default value 0.7 used." << endl;
			}
		}
		if (!strcmp(argv[i], "-lmd")) {
			float lmd = stof(argv[i + 1]);
			if (lmd > 0) {
				output.fista_lambda = lmd;
			}
			else {
				cout << "The lambda parameter of FISTA "
					"provided is in the wrong "
					"format (must be in (0, 1)). "
					"So, the default value 0.01 used." << endl;
			}
		}
		if (!strcmp(argv[i], "-ridg")) {
			output.output_ridgelets = argv[i + 1];
			out1 = true;
		}
		if (!strcmp(argv[i], "-sr")) {
			output.signal_recon = argv[i + 1];
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
		if (!strcmp(argv[i], "-nth")) {
			float n_threads = stof(argv[i + 1]);
			if (n_threads == floor(n_threads) && n_threads > 0) {
				output.nth = n_threads;
			}
			else {
				cout << "The value for number of threads "
					"provided is in the wrong format "
					"(must be a positive integer). "
					"So, will be computed automatically." << endl;
			}
		}
	}
	if (!inp1 || !out1) {
		cerr << "Please, provide at least one input AND one output file names" << endl;
		return EXIT_FAILURE;
	}
	return 0;
}

void DATA_SOURCE::short_summary(input_parse& params) {
	// Show summary on parameters will be used during computations

	cout << "Summary on parameters" << endl;
	cout << "-----------------------------------" << endl;
	cout << "Spherical ridgelets J: " << params.sph_J << endl;
	cout << "Spherical ridgelets rho: " << params.sph_rho << endl;
	cout << "Icosahedron tesselation order: " << params.lvl << endl;
	cout << "Maxima ODF threshold: " << params.max_odf_thresh << endl;
	cout << "FISTA lambda parameter: " << params.fista_lambda << endl;
	cout << "Number of splits: " << params.n_splits << endl;
	cout << "File(s) compression enabled: ";
	params.is_compress ? cout << "yes" : cout << "no";
	cout << endl << "-----------------------------------" << endl;
}

int DATA_SOURCE::readMask(string inputMask, MaskImagePointer& image) {
	// We need mask within the main program so it is implemented in that class

	// Temporary register factories cause don't use cmake
	itk::NrrdImageIOFactory::RegisterOneFactory();

	// Another way to store image data
	MaskReaderType::Pointer reader = MaskReaderType::New();

	// Make some inputfiles checks
	if (!is_path_exists(inputMask)) {
		cout << "Input mask image is not provided. Please, stop program and provide mask file "
			"if you forget to include it." << endl;
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
	// Commented header information might be utilized in future
	//dest->SetMetaDataDictionary(src->GetMetaDataDictionary());
	//dest->SetNumberOfComponentsPerPixel(h.comp_h);
	dest->SetSpacing(h.spc_h);
	dest->SetDirection(h.dirs_h);
	dest->SetOrigin(h.orig_h);
	dest->SetRegions(h.reg_h);
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
	SPH_RIDG<precisionType, MatrixType, VectorType> test;
	//cout << test.M0;
	UtilMath<precisionType, MatrixType> UM = UtilMath<precisionType, MatrixType>();
	MatrixType g;
	g = MatrixType::Zero(100, 3);
	UM.spiralsample(g, 2, 100);
	cout << g << endl;
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

void DATA_SOURCE::data_saving_info_out(unsigned long int coef_size, string name) {
	unsigned long long int orig_img_size = h.reg_h.GetSize()[0] * h.reg_h.GetSize()[1] * h.reg_h.GetSize()[2];
	unsigned long long int ridg_save = orig_img_size * coef_size * sizeof(precisionType);
	cout << "Also you need " << ridg_save / pow(1024, 3) << " GB of RAM to save " << name << endl;
}

void DATA_SOURCE::estimate_memory(MatrixType& s, MatrixType& A, input_parse& params) {
	unsigned n_splits = params.n_splits;

	// Get the number of threads
	unsigned n_threads = Eigen::nbThreads();
	cout << "The number of available threads: " << n_threads << endl;

	// Estimate dMRI Eigen matrix size
	unsigned long long int dmri_memory = s.size() * sizeof(precisionType);

	// Estimate memory consumption by FISTA solver
	unsigned long long int fista_memory_x = s.cols() * A.cols() * sizeof(precisionType);
	unsigned long long int fista_memory_loop = n_threads * 4 * (s.cols() / n_splits) * A.cols() * sizeof(precisionType);
	unsigned long long int total = dmri_memory + fista_memory_x + fista_memory_loop;

	cout << "IMPORTANT! To successfully finish computations you need approximately ";
	cout << total / pow(1024, 3) << " GB of RAM and virtual memory combined to compute spherical ridgelets." << endl;
	
	// Estimate memory consumption to save ridgelets
	if (!params.output_ridgelets.empty())
		data_saving_info_out(A.cols(), "ridgelets coefficients");

	// Estimate memory consumption to save reconstructed signal
	if (!params.signal_recon.empty())
		data_saving_info_out(A.rows(), "reconstructed signal");

	// Estimate memory consumption to save ODF max values and directions
	if (!params.output_fiber_max_odf.empty())
		data_saving_info_out(24, "ODF max");

	cout << "If you want to optimize memory consumption and computation speed, feel free to "
		"experiment with split coefficient (-nspl parameter). " << endl << endl;
}

int DATA_SOURCE::compute_splits(unsigned s_size) {
	unsigned n_threads = Eigen::nbThreads();
	int s = s_size / (n_threads * 2);
	cout << "An optimal number of splits for your dMRI image and CPU is " << s << endl;
	return s;
}

int DATA_SOURCE::DWI2Matrix(string &dmri_file, MaskImagePointer &mask, MatrixType &signal, MatrixType &grad_dirs)
{
	SIGNAL_GENERATOR sg(dmri_file);
	int res_dmri = sg.ExtractMatrix(mask, signal, grad_dirs);
	h = sg.h;
	if (res_dmri)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}

void DATA_SOURCE::Matrix2DWI(DiffusionImagePointer &img, MaskImagePointer &mask, MatrixType &arr) {
	int n_of_components = arr.rows();

	// Make a vector of zeros
	VariableVectorType zeros_vec;
	zeros_vec.SetSize(n_of_components);
#pragma omp parallel for
	for (int i = 0; i < n_of_components; ++i) {
		zeros_vec[i] = 0;
	}
	img->FillBuffer(zeros_vec);

	VariableVectorType vec_to_fill;
	vec_to_fill.SetSize(n_of_components);

	Iterator it(img, img->GetRequestedRegion());

	it.SetDirection(0);
	it.GoToBegin();
	unsigned vox = 0;

	if (mask) {
		// Mask iterator
		MaskIterator it_m(mask, mask->GetRequestedRegion());

		it_m.SetDirection(0);
		it_m.GoToBegin();

		while (!it.IsAtEnd())
		{
			while (!it.IsAtEndOfLine())
			{
				if (it_m.Get() == 1) {
					for (int i = 0; i < n_of_components; ++i) {
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
				for (int i = 0; i < n_of_components; ++i) {
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
	if (file.is_open()) {
		for (int i = 0; i < matrix.rows(); ++i) {
			for (int j = 0; j < matrix.cols(); ++j) {
				file << matrix(i, j) << " ";
			}
			file << endl;
		}
	}
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
		stringstream inner_stream(line);

		for (unsigned j = 0; j < number_of_cols; ++j)
			inner_stream >> matrix(i, j);
	}
}

template <typename T>
void DATA_SOURCE::printVec(const string& name, vector<T>& v) {
	cout << name << endl;
	for (auto it = v.cbegin(); it != v.cend(); ++it)
		std::cout << *it << ' ';

	cout << endl;
}

// Necessary templates for utilized types. Needed because implementation and declaration of function in separate files
template void DATA_SOURCE::printVec<int>(const string&, vector<int>&);
template void DATA_SOURCE::printVec<unsigned>(const string&, vector<unsigned>&);
template void DATA_SOURCE::printVec<Eigen::Index>(const string&, vector<Eigen::Index>&);

template int DATA_SOURCE::save_to_file<DiffusionImageType>(const string&, DiffusionImageType::Pointer&, bool);
