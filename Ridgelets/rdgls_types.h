#ifndef RDGLS_TYPES_H
#define RDGLS_TYPES_H

//libraries
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkNrrdImageIOFactory.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkImageLinearConstIteratorWithIndex.h>

#include <Eigen/Dense>

#include <omp.h>
#include <chrono>
#include <thread>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
/*
typedef itk::Vector<signed short, 59>   VectorType;
typedef itk::Image<VectorType, 3>	    DiffusionImageType;
typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
*/
typedef double                          VectorType;
typedef itk::VectorImage<VectorType, 3>	DiffusionImageType;
typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
typedef itk::ImageLinearConstIteratorWithIndex<DiffusionImageType> ConstIterator;

using namespace Eigen;
using namespace std;
using namespace std::chrono;

#endif
