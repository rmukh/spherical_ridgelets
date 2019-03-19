find_package(ITK REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(DEFINED Eigen3_DIR AND NOT EXISTS ${Eigen3_DIR})
  message(FATAL_ERROR "Eigen3_DIR variable is defined but corresponds to non-existing directory (${Eigen3_DIR})")
endif()

include(${ITK_USE_FILE})

include_directories(Ridgelets)
add_executable(sphridg 
Ridgelets/Ridgelets.cpp
Ridgelets/convhull_3d.h
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
target_link_libraries(sphridg ${ITK_LIBRARIES})
target_link_libraries(sphridg Eigen3::Eigen)
