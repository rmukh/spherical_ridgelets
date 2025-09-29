#---------------------------------------------------------------------------
# Get and build Eigen

set(EIGEN_TAG "3.4.1")

option(eigen_PATCH_FILE "Patch file to apply to eigen" OFF)
set(eigen_PATCH_COMMAND "")
if(eigen_PATCH_FILE)
    set(eigen_PATCH_COMMAND PATCH_COMMAND patch -p0 -i "${eigen_PATCH_FILE}")
endif(eigen_PATCH_FILE)

ExternalProject_Add(Eigen3
    PREFIX Eigen3
    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
    GIT_TAG "${EIGEN_TAG}"
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen/install
    ${eigen_PATCH_COMMAND}
    LOG_DOWNLOAD ON
    LOG_OUTPUT_ON_FAILURE ON
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE:STRING=Release
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
        -DCMAKE_INSTALL_PREFIX:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/eigen/install
        -DBUILD_TESTING:BOOL=OFF
        -DEIGEN_BUILD_PKGCONFIG:BOOL=OFF
        -DEIGEN_TEST_NOQT:BOOL=ON
)

set(Eigen3_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen/install)
