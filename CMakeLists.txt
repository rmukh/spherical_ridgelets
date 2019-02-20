cmake_minimum_required (VERSION 3.5)

#use c++ 11
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project (Ridgelets)

# Command to output information to the console (debugging)
message ("cxx Flags:" ${CMAKE_CXX_FLAGS})

add_executable(sphridg Ridgelets/Ridgelets.cpp)