cmake_minimum_required (VERSION 3.13)

#use c++ 11
set(_msg "Setting C++ standard")
message(STATUS "${_msg}")
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")

# Command to output information to the console (debugging)
message ("cxx Flags:" ${CMAKE_CXX_FLAGS})

if(NOT DEFINED ALL_LIBS_SET)
    set(ALL_LIBS_SET false)
endif()

set(PRJ_NAME Ridgelets)

if(NOT ALL_LIBS_SET)
    message("Start buildining and configuring all dependencies")
    project(SuperBuild_ALL_LIBS)
    include("${CMAKE_CURRENT_SOURCE_DIR}/libs.cmake")
    return()
else()
    project (${PRJ_NAME})
    include("${CMAKE_CURRENT_SOURCE_DIR}/${PRJ_NAME}.cmake")
    return()
endif()
