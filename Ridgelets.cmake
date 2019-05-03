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

include_directories(Ridgelets)
add_executable(sphridg 
Ridgelets/Ridgelets.cpp
Ridgelets/convhull_3d.h
Ridgelets/SIGNAL_GENERATOR.cpp
Ridgelets/SIGNAL_GENERATOR.h
Ridgelets/DATA_SOURCE.cpp
Ridgelets/DATA_SOURCE.h
Ridgelets/packages.config
Ridgelets/rdgls_types.h
Ridgelets/Ridgelets.cpp
Ridgelets/SOLVERS.cpp
Ridgelets/SOLVERS.h
Ridgelets/SPH_RIDG.cpp
Ridgelets/SPH_RIDG.h
Ridgelets/UtilMath.cpp
Ridgelets/UtilMath.h)
target_link_libraries(sphridg Eigen3::Eigen)
target_link_libraries(sphridg ${ITK_LIBRARIES})
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
