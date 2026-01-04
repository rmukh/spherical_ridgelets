# ExternalProject to fetch cxxopts (header-only)
# Usage: include(this file) when BUILD_CLI is ON.
include(ExternalProject)

set(CXXOPTS_TAG "v3.0.0" CACHE STRING "Tag/branch of cxxopts to fetch")
ExternalProject_Add(cxxopts
  GIT_REPOSITORY "https://github.com/jarro2783/cxxopts"
  GIT_TAG "${CXXOPTS_TAG}"
  SOURCE_DIR "${CMAKE_BINARY_DIR}/cxxopts-src"
  BINARY_DIR "${CMAKE_BINARY_DIR}/cxxopts-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  LOG_OUTPUT_ON_FAILURE ON
)

# Make the header available via ${CMAKE_BINARY_DIR}/cxxopts-src/include
# Consumers should add that directory to include paths:
# target_include_directories(<target> PRIVATE ${CMAKE_BINARY_DIR}/cxxopts-src/include)
