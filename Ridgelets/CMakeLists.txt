cmake_minimum_required(VERSION 3.6)
include("${CMAKE_CURRENT_SOURCE_DIR}/../Common.cmake")

if(DEFINED BUILD_SHARED_LIBS)
    message(STATUS "Building the library of type: " ${BUILD_SHARED_LIBS})
    set(LIB_BUILD_TYPE ${BUILD_SHARED_LIBS})
else()
    message(STATUS "Building the library with default type: STATIC")
    set(LIB_BUILD_TYPE STATIC)
endif()

project(${lib_name})

add_library(
${lib_name} ${LIB_BUILD_TYPE}
SIGNAL_GENERATOR.cpp
DATA_SOURCE.cpp
packages.config
SOLVERS.cpp
SPH_RIDG.cpp
UtilMath.cpp
)

set_target_properties(
${lib_name} PROPERTIES PUBLIC_HEADER 
Ridgelets/convhull_3d.h
Ridgelets/SIGNAL_GENERATOR.h
Ridgelets/DATA_SOURCE.h
Ridgelets/rdgls_types.h
Ridgelets/SOLVERS.h
Ridgelets/SPH_RIDG.h
Ridgelets/UtilMath.h
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ITK REQUIRED)
find_package(OpenMP)

if(DEFINED Eigen3_DIR AND NOT EXISTS ${Eigen3_DIR})
  message(FATAL_ERROR "Eigen3_DIR variable is defined but corresponds to non-existing directory (${Eigen3_DIR})")
endif()

if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory (${ITK_DIR})")
endif()

include(${ITK_USE_FILE})

target_link_libraries(${lib_name} Eigen3::Eigen ${ITK_LIBRARIES}) 

add_executable(sphridg Ridgelets.cpp)
target_link_libraries(sphridg Eigen3::Eigen ${ITK_LIBRARIES} ${lib_name})

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP will be linked")
  target_link_libraries(sphridg OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP can't be linked. Probably you use CMake < 3.9")
  message(STATUS "Attempt to link Threads package...")
endif()

if(NOT TARGET OpenMP::OpenMP_CXX)
  find_package(Threads REQUIRED)
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX
               PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
  # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
  set_property(TARGET OpenMP::OpenMP_CXX
               PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()
target_link_libraries(sphridg OpenMP::OpenMP_CXX)