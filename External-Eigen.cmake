#---------------------------------------------------------------------------
# Get and build Eigen

set(EIGEN_TAG "for/itk-20231108-master-4d54c43d")
set(EIGEN3_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen/install)
set(EIGEN3_CMAKE_DIR ${EIGEN3_INSTALL_DIR}/share)
option(eigen_PATCH_FILE "Patch file to apply to eigen" OFF)
set(eigen_PATCH_COMMAND "")
if(eigen_PATCH_FILE)
    set(eigen_PATCH_COMMAND PATCH_COMMAND patch -p0 -i "${eigen_PATCH_FILE}")
endif(eigen_PATCH_FILE)

ExternalProject_Add(Eigen3
    PREFIX Eigen3
    GIT_REPOSITORY "https://github.com/InsightSoftwareConsortium/eigen"
    GIT_TAG "${EIGEN_TAG}"
    UPDATE_DISCONNECTED ON
    INSTALL_DIR ${EIGEN3_INSTALL_DIR}
    ${eigen_PATCH_COMMAND}
    LOG_DOWNLOAD ON
    LOG_OUTPUT_ON_FAILURE ON
    CMAKE_ARGS
        ${SPH_EXTERNAL_CMAKE_ARGS}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DEIGEN_BUILD_TESTING:BOOL=OFF
        -DEIGEN_BUILD_PKGCONFIG:BOOL=OFF
        -DEIGEN_BUILD_DOC:BOOL=OFF
        -DEIGEN_TEST_NOQT:BOOL=ON
    CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${EIGEN3_INSTALL_DIR}
)

set(Eigen3_DIR ${EIGEN3_CMAKE_DIR})
