# Sanity checks
find_package(ITK REQUIRED)

if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory")
endif()

include(${ITK_USE_FILE})

add_executable(sphridg Ridgelets/Ridgelets.cpp)
target_link_libraries(sphridg ${ITK_LIBRARIES})
