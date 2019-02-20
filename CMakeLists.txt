cmake_minimum_required (VERSION 2.6)
 
project (Ridgelets)

# Command to output information to the console
# Useful for displaying errors, warnings, and debugging
message ("cxx Flags:" ${CMAKE_CXX_FLAGS})
 
add_executable(Ridgelets Ridgelets.cpp)