#include "pch.h"
#include "DATA_SOURCE.h"

int DATA_SOURCE::readNRRD(
	string inputVolume, MatrixType &GradientDirections,
	DiffusionImagePointer &image, unsigned &nGradImgs,
	unsigned &nOfImgs) {
	bool is_b0 = false;
	double b0 = 0;
	double x, y, z;

	//temporary register factories cause don't use cmake
	itk::NrrdImageIOFactory::RegisterOneFactory();

	//Another way to store image data

	typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	//get image

	reader->SetFileName(inputVolume);
	try
	{
		reader->Update();
		image = reader->GetOutput();
	}
	catch (itk::ExceptionObject & err)
	{
		cerr << "ExceptionObject caught !" << endl;
		cerr << err << endl;

		// Since the goal of the example is to catch the exception,
		// we declare this a success.
		return EXIT_SUCCESS;
	}

	//Get and process header information
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
			sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);

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
	//Normalize directions
	GradientDirections.rowwise().normalize();

	cout << "Number of gradient images: " << nGradImgs << " and Number of reference images: " << nOfImgs - nGradImgs << endl;
	cout << "b value " << b0 << endl;
	if (!is_b0)
	{
		cerr << "BValue not specified in header file" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
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
	MatrixType A = test.RBasis(g);
	MatrixType normA = test.normBasis(A);

	test.normBasis(normA);
}

MatrixType DATA_SOURCE::readF(string f, int x, int y) {
	ifstream file(f);
	string str;
	MatrixType v(x, y);
	int g = 0;
	while (getline(file, str))
	{
		vector<string> result;
		istringstream iss(str);
		int i = 0;
		for (string s; iss >> s;) {
			v(g, i) = atof(s.c_str());
			i++;
		}
		g++;
	}
	return v;
}

void DATA_SOURCE::readTestData(MatrixType& g, MatrixType& s) {
	g = readF("D:\\test\\bvec.txt", 51, 3);
	s = readF("D:\\test\\sign.txt", 51, 16);

	//normalize
	MatrixType gnorm = g.array().pow(2).rowwise().sum().sqrt();
	for (int i = 0; i < 51; ++i)
		g.row(i) = g.row(i) / gnorm(i);
}

void DATA_SOURCE::DWI2Matrix(DiffusionImagePointer &img, MatrixType &signal, unsigned &nGradImgs, unsigned &nOfImgs) {
	DiffusionImageType::SizeType sz = img->GetLargestPossibleRegion().GetSize();
	unsigned first_grad_image_index = nOfImgs - nGradImgs;
	int N_of_voxels = sz[0] * sz[1] * sz[2];
	signal = MatrixType::Zero(nGradImgs, N_of_voxels);
	DiffusionImageType::PixelType voxel_content;

	//Iterate over all voxels
	ConstIterator it(img, img->GetRequestedRegion());
	unsigned vox = 0;
	it.SetDirection(0);
	it.GoToBegin();
	while (!it.IsAtEnd())
	{
		while (!it.IsAtEndOfLine())
		{
			voxel_content = it.Get();
			for (unsigned i = first_grad_image_index; i < nOfImgs; ++i)
				signal(i - first_grad_image_index, vox) = voxel_content.GetElement(i);

			++vox;
			++it;
		}
		it.NextLine();
	}
}

void DATA_SOURCE::matrixToFile(const string& fname, MatrixType& matrix) {
	ofstream file(fname);
	if (file.is_open())
		file << matrix << '\n';
}

void DATA_SOURCE::fileToMatrix(const string& fname, MatrixType& matrix)
{
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

template void DATA_SOURCE::printVec<int>(const string&, vector<int>&); 
template void DATA_SOURCE::printVec<Eigen::Index>(const string&, vector<Eigen::Index>&);