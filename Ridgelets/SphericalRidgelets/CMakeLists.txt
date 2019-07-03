cmake_minimum_required(VERSION 3.6)
if(BUILD_TU)
	include("${SPH_DIR}/../../Common.cmake")
endif()
project(SphericalRidgelets VERSION 1.2)
include(GNUInstallDirs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ITK REQUIRED)
find_package(OpenMP)

if(DEFINED BUILD_SHARED_LIBS)
    message(STATUS "Building the library of type: " ${BUILD_SHARED_LIBS})
    set(LIB_BUILD_TYPE ${BUILD_SHARED_LIBS})
else()
    message(STATUS "Building the library with default type: STATIC")
    set(LIB_BUILD_TYPE STATIC)
endif()

add_library(
Spherical_Ridgelets ${LIB_BUILD_TYPE}
${SPH_DIR}/src/SIGNAL_GENERATOR.cpp
${SPH_DIR}/src/DATA_SOURCE.cpp
${SPH_DIR}/src/convhull_3d.cpp
${SPH_DIR}/include/convhull_3d.h
${SPH_DIR}/include/UtilMath.h
${SPH_DIR}/include/SIGNAL_GENERATOR.h
${SPH_DIR}/include/DATA_SOURCE.h
${SPH_DIR}/include/rdgls_types.h
${SPH_DIR}/include/SOLVERS.h
${SPH_DIR}/include/SPH_RIDG.h
${SPH_DIR}/include/SOLVERS.hpp
${SPH_DIR}/include/SPH_RIDG.hpp
${SPH_DIR}/include/UtilMath.hpp
)

set_target_properties(Spherical_Ridgelets
    PROPERTIES
      OUTPUT_NAME "Spherical_Ridgelets"
)

add_library(Spherical_Ridgelets::Spherical_Ridgelets ALIAS Spherical_Ridgelets)

#target_compile_features(Spherical_Ridgelets
#    PUBLIC
#      cxx_attributes
#      cxx_defaulted_functions
#      cxx_deleted_functions
#      cxx_final
#)

target_include_directories(Spherical_Ridgelets
    PUBLIC
    	"$<BUILD_INTERFACE:${SPH_DIR}/src;${SPH_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

include(${ITK_USE_FILE})
target_link_libraries(Spherical_Ridgelets Eigen3::Eigen)
target_link_libraries(Spherical_Ridgelets ${ITK_LIBRARIES})

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP will be linked")
  target_link_libraries(Spherical_Ridgelets OpenMP::OpenMP_CXX)
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
target_link_libraries(Spherical_Ridgelets OpenMP::OpenMP_CXX)

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SphericalRidgelets)

add_executable(sphridg ${SPH_DIR}/../Ridgelets.cpp)
target_link_libraries(sphridg PUBLIC Eigen3::Eigen ${ITK_LIBRARIES} OpenMP::OpenMP_CXX Spherical_Ridgelets)

include(CMakePackageConfigHelpers)

install(
    TARGETS
      Spherical_Ridgelets
    EXPORT
      SphericalRidgeletsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT SphericalRidgeletsTargets
    FILE SphericalRidgeletsTargets.cmake
    NAMESPACE Spherical_Ridgelets::
    DESTINATION ${INSTALL_CONFIGDIR}
)

install(DIRECTORY ${SPH_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h*"
)

# ConfigVersion file

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SphericalRidgeletsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/SphericalRidgeletsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SphericalRidgeletsConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install all the helper files
install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/SphericalRidgeletsConfig.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/SphericalRidgeletsConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
