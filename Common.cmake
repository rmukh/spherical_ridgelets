# Use c++ 11
set(_msg "Setting C++ standard")
message(STATUS "${_msg}")
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    if(MSVC)
        string(APPEND CMAKE_CXX_FLAGS " /Zc:__cplusplus")
    endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")

# Enable optimization
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES RELEASE)
    if (MSVC)
        include(CheckCXXCompilerFlag)

        add_compile_options(/O2)

        check_cxx_compiler_flag("/arch:AVX2" COMPILER_OPT_ARCH_AVX2_SUPPORTED)
        if(COMPILER_OPT_ARCH_AVX2_SUPPORTED)
          add_compile_options(/arch:AVX2)
        else()
          check_cxx_compiler_flag("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
          if(COMPILER_OPT_ARCH_AVX_SUPPORTED)
           add_compile_options(/arch:AVX)
          endif()
          check_cxx_compiler_flag("/arch:SSE2" COMPILER_OPT_SSE2_SUPPORTED)
          if(COMPILER_OPT_SSE2_SUPPORTED)
           add_compile_options(/arch:SSE2)
          endif()
        endif()

        if (JUST_BUILD)
            add_compile_options(/GL /GA)
        endif()
    else()
        add_compile_options(-O3)
        add_compile_options(-march=native)
    endif()

    message(STATUS "O3 level optimization and native CPU features support enabled")
endif()

message("Build type: ${CMAKE_BUILD_TYPE}")

# Check if correct build type passed
if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
    message(FATAL_ERROR "Invalid value for CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# Command to output information to the console (debugging)
message("cxx Flags:" ${CMAKE_CXX_FLAGS})
