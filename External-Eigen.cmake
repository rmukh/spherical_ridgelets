#---------------------------------------------------------------------------
# Get and build Eigen

set(EIGEN_TAG "3.3.7")
ExternalProject_Add(Eigen3
  GIT_REPOSITORY "${git_protocol}://github.com/eigenteam/eigen-git-mirror"
  GIT_TAG "${EIGEN_TAG}"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen
  INSTALL_COMMAND ""
)

set(Eigen3_DIR_INSTALLED ${CMAKE_BINARY_DIR}/eigen/src/Eigen3/cmake)

