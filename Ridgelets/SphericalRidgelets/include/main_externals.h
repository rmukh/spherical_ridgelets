#ifndef MAIN_EXTERNSLS_H
#define MAIN_EXTERNSLS_H

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

#ifdef _OPENMP
  #include <omp.h>
#endif
#include <chrono>
#include <thread>
#include <unordered_map>
#include <sys/stat.h>
#include <numeric>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

// Necessary namespaces
using namespace Eigen;
using namespace std;
using namespace std::chrono;

#endif
