// Ridgelets.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//precompiled header
#include "pch.h"

//user defined classes
#include "UtilMath.h"
#include "SOLVERS.h"
#include "SPH_RIDG.h"
#include "DATA_SOURCE.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkPointSource.h>

int main()
{
	DATA_SOURCE data;
	/*
	MatrixType GradientDirections(0, 3);
	DiffusionImagePointer image;
	unsigned nGradImgs = 0;
	unsigned nOfImgs = 0;
	data.readNRRD("C:\\Users\\renat\\Desktop\\01009-dwi-Ed.nhdr", GradientDirections, image, nGradImgs, nOfImgs);

	const DiffusionImageType::IndexType pixelIndex = { {27,29,37} }; //Position {X,Y,Z}
	DiffusionImageType::PixelType value = image->GetPixel(pixelIndex);

	cout << "Gradient Directions: \n" << GradientDirections << endl;
	cout << value << endl;

	//4D dMRI image to Eigen 2D Matrix
	MatrixType signal;
	data.DWI2Matrix(image, signal, nGradImgs, nOfImgs);
	image = nullptr;

	//cout << "voxel 50000" << signal.col(50000) << endl;
	*/
	MatrixType g(51, 3);
	MatrixType s(51, 16);

	data.readTestData(g, s);

	//icosahedron
	UtilMath m;
	int level = 4;

	double C = 1 / sqrt(1.25);
	MatrixType t = (2 * m.PI / 5.0) * VectorXd::LinSpaced(5, 0, 4);
	MatrixType u1(5, 3);
	u1 << C * t.array().cos(), C * t.array().sin(), C * 0.5 * MatrixType::Ones(5, 1);
	MatrixType u2(5, 3);
	u2 << C * (t.array() + 0.2 * m.PI).cos(), C * (t.array() + 0.2 * m.PI).sin(), -0.5 * C * MatrixType::Ones(5, 1);
	MatrixType u(12, 3);
	u << 0, 0, 1, u1, u2, 0, 0, -1;

	/*
	   Convex hull part (might be conflict between VTK 7 and 8 versions)
	*/

	// Convert from Eigen MatrixType to vtkPoints (probably make as an function)
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (unsigned i = 0; i < u.rows(); ++i)
		points->InsertNextPoint(u(i, 0), u(i, 1), u(i, 2));

	// Dataset to represent verticies (points)
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);

	// Create the convex hull of the pointcloud
	vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay->SetInputData(polydata);
	delaunay->Update();

	// Convert to polygonal type
	vtkSmartPointer<vtkUnstructuredGrid> raw = delaunay->GetOutput();
	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	geometryFilter->SetInputData(raw);
	geometryFilter->Update();

	vtkSmartPointer<vtkPolyData> polyPoints = geometryFilter->GetOutput();

	unsigned NCells = polyPoints->GetNumberOfCells();
	cout << "Output has " << NCells << " cells." << endl;

	// Convert to Eigen matrix
	MatrixType fcs(NCells, 3);
	for (vtkIdType i = 0; i < NCells; ++i) {
		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		polyPoints->GetCellPoints(i, cellPointIds);

		for (vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); ++j)
			fcs(i, j) = cellPointIds->GetId(j);
	}

	cout << "faces" << endl << fcs.array() + 1;
	//for (vtkIdType i = 0; i < raw->GetNumberOfPoints(); i++)
	//{
	//	double p[3];
	//	raw->GetCell(i, p);
	//	cout << "point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << endl;
	//}

	//
	//high_resolution_clock::time_point t1 = high_resolution_clock::now(); //start timer point

	//SPH_RIDG ridg(2, 0.5);
	//MatrixType A = ridg.RBasis(g);
	//A = ridg.normBasis(A)

	//SOLVERS slv(A, s, 0.1);
	//MatrixType C = slv.FISTA();
	//cout << endl << C;

	////ODF
	////MatrixType Q = ridg.QBasis(nu); //Build a Q basis
	////MatrixType ODF = C * Q.transpose(); //Computer ODF

	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//auto duration = duration_cast<seconds>(t2 - t1).count();
	//cout << "Execution time " << duration << " seconds" << endl;

	return 0;
}
