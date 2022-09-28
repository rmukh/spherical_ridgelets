#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "SPH_RIDG.h"

namespace py = pybind11;
using namespace pybind11::literals;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorType;


PYBIND11_MODULE(Spherical_Ridgelets_py, m)
{
    m.doc() = R"pbdoc(
        Spherical Ridgelets
        -----------------------

        .. currentmodule:: Spherical_Ridgelets_py

        .. autosummary::
           :toctree: _generate

           SphericalRidgelets
    )pbdoc";
    py::class_<SPH_RIDG<double, MatrixType, VectorType>>(m, "SPH_RIDG");
}