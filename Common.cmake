# Use c++ 11
set(_msg "Setting C++ standard")
message(STATUS "${_msg}")
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")

# Enable optimization
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES RELEASE)
    add_compile_options(-O3)
    add_compile_options(-march=native)
    message(STATUS "O3 level optimization and native CPU features support enabled")
endif()

# Check if correct build type passed
if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
    message(FATAL_ERROR "Invalid value for CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# Command to output information to the console (debugging)
message("cxx Flags:" ${CMAKE_CXX_FLAGS})
