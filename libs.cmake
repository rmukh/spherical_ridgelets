set(_SPH_DOWNLOAD_EIGEN OFF)
set(_SPH_DOWNLOAD_ITK OFF)

if(DEFINED Eigen3_DIR AND NOT "${Eigen3_DIR}" STREQUAL "")
  file(TO_CMAKE_PATH "${Eigen3_DIR}" Eigen3_DIR)
  if(NOT EXISTS "${Eigen3_DIR}/Eigen3Config.cmake" AND
     NOT EXISTS "${Eigen3_DIR}/eigen3-config.cmake")
    message(FATAL_ERROR
      "Eigen3_DIR does not contain Eigen3Config.cmake: ${Eigen3_DIR}\n"
      "Leave Eigen3_DIR empty to download the pinned Eigen version.")
  endif()
else()
  unset(Eigen3_DIR CACHE)
  unset(Eigen3_DIR)
  set(_SPH_DOWNLOAD_EIGEN ON)
endif()

if(DEFINED ITK_DIR AND NOT "${ITK_DIR}" STREQUAL "")
  file(TO_CMAKE_PATH "${ITK_DIR}" ITK_DIR)
  if(NOT EXISTS "${ITK_DIR}/ITKConfig.cmake")
    message(FATAL_ERROR
      "ITK_DIR does not contain ITKConfig.cmake: ${ITK_DIR}\n"
      "Leave ITK_DIR empty to download the pinned ITK version.")
  endif()
  file(STRINGS "${ITK_DIR}/ITKConfig.cmake" _SPH_ITK_VERSION_MAJOR_LINE
    REGEX "set\\(ITK_VERSION_MAJOR")
  string(REGEX REPLACE ".*ITK_VERSION_MAJOR \"?([0-9]+)\"?.*" "\\1"
    _SPH_ITK_VERSION_MAJOR "${_SPH_ITK_VERSION_MAJOR_LINE}")
  if(_SPH_ITK_VERSION_MAJOR AND _SPH_ITK_VERSION_MAJOR VERSION_LESS 5)
    message(WARNING
      "Ignoring ITK_DIR because it points to ITK ${_SPH_ITK_VERSION_MAJOR}.x: ${ITK_DIR}\n"
      "This project requires ITK 5 or newer with the current MinGW/GCC toolchain. "
      "The pinned ITK version will be downloaded instead.")
    unset(ITK_DIR CACHE)
    unset(ITK_DIR)
    set(_SPH_DOWNLOAD_ITK ON)
  endif()
else()
  unset(ITK_DIR CACHE)
  unset(ITK_DIR)
  set(_SPH_DOWNLOAD_ITK ON)
endif()

if(_SPH_DOWNLOAD_EIGEN OR _SPH_DOWNLOAD_ITK)
  find_package(Git REQUIRED)
endif()

include(ExternalProject)

if(_SPH_DOWNLOAD_EIGEN)
  include(${CMAKE_SOURCE_DIR}/External-Eigen.cmake)
  message(STATUS "Using External Project for Eigen")
else()
  message(STATUS "Using specified Eigen3_DIR: ${Eigen3_DIR}")
  if(NOT TARGET Eigen3)
    add_custom_target(Eigen3)
  endif()
endif()

if(_SPH_DOWNLOAD_ITK)
  include(${CMAKE_SOURCE_DIR}/External-ITK.cmake)
  message(STATUS "Using External Project for ITK")
else()
  message(STATUS "Using specified ITK_DIR: ${ITK_DIR}")
  if(NOT TARGET ITK)
    add_custom_target(ITK)
  endif()
endif()
