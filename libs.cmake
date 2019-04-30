#Check if git installed
find_package(Git)
if(NOT GIT_FOUND)
  message(ERROR "Cannot find git. git is required for Superbuild")
endif()

option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

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

set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
if(APPLE)
  list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
  )
endif(APPLE)

include(ExternalProject)
# Find Eigen and build if necessary
if(Eigen3_DIR)
    set(Eigen3_DIR_INSTALLED ${Eigen3_DIR})
    message(STATUS "Using specified Eigen3_DIR: ${Eigen3_DIR}")
else(Eigen3_DIR)
    include(${CMAKE_SOURCE_DIR}/External-Eigen.cmake)
    message(STATUS "Using External Project for Eigen")
endif(Eigen3_DIR)

include(ExternalProject)
# Find ITK and build if necessary
if(ITK_DIR)
    set(ITK_DIR_INSTALLED ${ITK_DIR})
    message(STATUS "Using specified ITK_DIR: ${ITK_DIR}")
else(ITK_DIR)
    include(${CMAKE_SOURCE_DIR}/External-ITK.cmake)
    message(STATUS "Using External Project for ITK")
endif(ITK_DIR)

set(${PRJ_NAME}_DEPENDENCIES
  Eigen3
  ITK
  )

set(proj ${PRJ_NAME})
# Configure and build Ridgelets
ExternalProject_Add(${proj}
  DEPENDS ${${PRJ_NAME}_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${PRJ_NAME}-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DALL_LIBS_SET:BOOL=true
    -DEigen3_DIR:PATH=${Eigen3_DIR_INSTALLED}
    -DITK_DIR:PATH=${ITK_DIR_INSTALLED}
  INSTALL_COMMAND ""
  )

ExternalProject_Add_Step(${proj} forcebuild
    COMMAND ${CMAKE_COMMAND} -E remove
    ${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
    COMMENT "Forcing build step for '${proj}'"
    DEPENDEES build
    ALWAYS 1
)

