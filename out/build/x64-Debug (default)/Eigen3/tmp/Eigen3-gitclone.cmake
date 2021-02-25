
if(NOT "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitinfo.txt" IS_NEWER_THAN "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "C:/Program Files/Git/cmd/git.exe"  clone --no-checkout "git://github.com/eigenteam/eigen-git-mirror" "Eigen3"
    WORKING_DIRECTORY "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'git://github.com/eigenteam/eigen-git-mirror'")
endif()

execute_process(
  COMMAND "C:/Program Files/Git/cmd/git.exe"  checkout 3.3.7 --
  WORKING_DIRECTORY "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '3.3.7'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "C:/Program Files/Git/cmd/git.exe"  submodule update --recursive --init 
    WORKING_DIRECTORY "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitinfo.txt"
    "D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'D:/VSProjects/SphericalRidgelets/out/build/x64-Debug (default)/Eigen3/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt'")
endif()

