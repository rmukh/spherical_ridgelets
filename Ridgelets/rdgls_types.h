#ifndef RDGLS_TYPES_H
#define RDGLS_TYPES_H

// Libraries
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkPointSource.h>

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

// Necessary types defenitions
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
/* Alternate way to define dMRI volume
typedef itk::Vector<signed short, 59>   VectorType;
typedef itk::Image<VectorType, 3>	    DiffusionImageType;
typedef DiffusionImageType::Pointer	    DiffusionImagePointer;
*/
typedef double                          VectorType;
typedef itk::VectorImage<VectorType, 3>	DiffusionImageType;
typedef itk::VariableLengthVector<VectorType> VariableVectorType;
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

#endif
