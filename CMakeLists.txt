cmake_minimum_required(VERSION 3.22)
project(StandaloneSphericalRidgelets)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#Configure external project with the same cmake generator
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}" )
endif()

if(UKF_USE_FLOAT)
  set(USE_FLOAT 1)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/libs.cmake")

message("Build type: ${CMAKE_BUILD_TYPE}")
message("C options: ${CMAKE_C_FLAGS}")
message("Cxx options: ${CMAKE_CXX_FLAGS}")

if(JUST_BUILD)
  message(STATUS "Build Spherical Ridgelets as a standalone application")
  set(SPH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Ridgelets/SphericalRidgelets)
  set(PRJ_NAME SphericalRidgeletsStandalone)
  set(${PRJ_NAME}_DEPENDENCIES 
    Eigen3
    ITK
  )

  ExternalProject_Add(${PRJ_NAME}
    DEPENDS ${${PRJ_NAME}_DEPENDENCIES}
    DOWNLOAD_COMMAND ""
    INSTALL_COMMAND ""
    SOURCE_DIR ${SPH_DIR}
    BINARY_DIR ${PRJ_NAME}-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DEigen3_DIR:PATH=${Eigen3_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DSPH_DIR:PATH=${SPH_DIR}
      -DBUILD_TU:STRING=1
      -DUSE_FLOAT:STRING=${USE_FLOAT}
    )
else()
  message(STATUS "Include Spherical Ridgelets for ExternalProject_Add() in ${SPH_DIR}")
  include(${CMAKE_CURRENT_SOURCE_DIR}/Ridgelets/SphericalRidgelets/CMakeLists.txt)
endif()
