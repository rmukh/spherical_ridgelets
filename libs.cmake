if(NOT Eigen3_DIR OR NOT ITK_DIR)
    #Check if git installed
    find_package(Git)
    if(NOT GIT_FOUND)
        message(ERROR "Cannot find git. git is required to build project in your specific case")
    endif()

    option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

    set(git_protocol "git")
    if(NOT USE_GIT_PROTOCOL)
        set(git_protocol "http")
    endif()
endif()

include(ExternalProject)
# Find Eigen and build if necessary
if(Eigen3_DIR)
    message(STATUS "Using specified Eigen3_DIR: ${Eigen3_DIR}")
else(Eigen3_DIR)
    include(${CMAKE_SOURCE_DIR}/External-Eigen.cmake)
    message(STATUS "Using External Project for Eigen")
endif(Eigen3_DIR)

include(ExternalProject)
# Find ITK and build if necessary
if(ITK_DIR)
    message(STATUS "Using specified ITK_DIR: ${ITK_DIR}")
else(ITK_DIR)
    include(${CMAKE_SOURCE_DIR}/External-ITK.cmake)
    message(STATUS "Using External Project for ITK")
endif(ITK_DIR)