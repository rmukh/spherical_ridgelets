#ifndef RDGLS_TYPES_H
#define RDGLS_TYPES_H

// Libraries
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <itkNrrdImageIOFactory.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <Eigen/Dense>

#include <omp.h>
#include <chrono>
#include <thread>
#include <unordered_map>
#include <sys/stat.h>
#include <numeric>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

//#define USE_FLOAT 1
#if USE_FLOAT
typedef float precisionType;
#define CONVHULL_3D_USE_FLOAT_PRECISION 1
#else
typedef double precisionType;
#endif

// Necessary types defenitions
typedef Eigen::Matrix<precisionType, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
typedef Eigen::Matrix<precisionType, Eigen::Dynamic, 1> VectorType;
/* Alternate way to define dMRI volume
typedef itk::Vector<signed short, 59>   VectorType;
typedef itk::Image<VectorType, 3>	    DiffusionImageType;
typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
*/
typedef itk::VectorImage<precisionType, 3>	DiffusionImageType;
typedef itk::VariableLengthVector<precisionType> VariableVectorType;
typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
typedef itk::ImageFileReader<DiffusionImageType> ImageReaderType;
typedef itk::ImageLinearConstIteratorWithIndex<DiffusionImageType> ConstIterator;
typedef itk::ImageLinearIteratorWithIndex<DiffusionImageType> Iterator;

// For mask file
typedef unsigned char                   MaskPixelType;
typedef itk::Image<MaskPixelType, 3>    MaskImageType;
typedef MaskImageType::Pointer	        MaskImagePointer;
typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
typedef itk::ImageLinearConstIteratorWithIndex<MaskImageType> MaskIterator;

// Necessary namespaces
using namespace Eigen;
using namespace std;
using namespace std::chrono;

// Useful structures
struct dMRI_h_info
{
	DiffusionImageType::SpacingType spc_h;
	DiffusionImageType::DirectionType dirs_h;
	DiffusionImageType::PointType orig_h;
	DiffusionImageType::RegionType reg_h;
	unsigned comp_h;
};

#endif
