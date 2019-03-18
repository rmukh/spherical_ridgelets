#---------------------------------------------------------------------------
# Get and build itk

set( ITK_TAG "v4.13.1" )
ExternalProject_Add(ITK
  GIT_REPOSITORY "${git_protocol}://itk.org/ITK.git"
  GIT_TAG "${ITK_TAG}"
  SOURCE_DIR ITK
  BINARY_DIR ITK-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DBUILD_TESTING:BOOL=OFF
    -DBUILD_EXAMPLES:BOOL=OFF
    -DITK_LEGACY_REMOVE:BOOL=ON
    -DITKV3_COMPATIBILITY:BOOL=OFF
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF
    -DModule_ITKCommon:BOOL=ON
    -DModule_ITKIONIFTI:BOOL=ON
    -DModule_ITKIONRRD:BOOL=ON
    -DModule_MGHIO:BOOL=ON
    -DModule_ITKIOMINC:BOOL=ON
    -DModule_ITKIOXML:BOOL=ON # For SlicerExecutionModel
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
    -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
    -DITK_WRAPPING:BOOL=OFF #${BUILD_SHARED_LIBS} ## HACK:  QUICK CHANGE
    -DITK_USE_FFTWD:BOOL=OFF
    -DITK_USE_FFTWF:BOOL=OFF
  INSTALL_COMMAND ""
)

set(ITK_DIR ${CMAKE_BINARY_DIR}/ITK-build)
