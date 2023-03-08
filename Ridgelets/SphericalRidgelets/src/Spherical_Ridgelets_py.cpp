#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "SPH_RIDG.h"
#include "SOLVERS.h"

namespace py = pybind11;
using namespace pybind11::literals;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorType;


PYBIND11_MODULE(Spherical_Ridgelets_py, m)
{
    m.doc() = R"pbdoc(
        Spherical Ridgelets binding
        -----------------------

        .. currentmodule:: Spherical_Ridgelets_py

        .. autosummary::
           :toctree: _generate

           SphericalRidgelets
    )pbdoc";
    py::class_<SPH_RIDG<double, MatrixType, VectorType>>(m, "SPH_RIDG")
        .def(py::init())
        .def(py::init<unsigned, double>())
        .def("init", &SPH_RIDG<double, MatrixType, VectorType>::init)
        .def("GetRBasis", &SPH_RIDG<double, MatrixType, VectorType>::GetRBasis_py, py::return_value_policy::reference_internal)
        .def("normBasis", &SPH_RIDG<double, MatrixType, VectorType>::normBasis, py::return_value_policy::reference_internal)
        .def("QBasis", &SPH_RIDG<double, MatrixType, VectorType>::QBasis, py::return_value_policy::reference_internal);
    
    py::class_<SOLVERS<double, MatrixType, MatrixType>>(m, "SOLVERS")
        .def(py::init())
        .def(py::init<MatrixType&, MatrixType&>())
        .def(py::init<MatrixType&, MatrixType&, double>())
        .def("SolveSingleSignal", &SOLVERS<double, MatrixType, MatrixType>::SolveSingle_py);
}