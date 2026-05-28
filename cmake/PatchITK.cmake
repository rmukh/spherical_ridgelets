if(NOT DEFINED ITK_SOURCE_DIR)
  message(FATAL_ERROR "ITK_SOURCE_DIR is required")
endif()

set(_header "${ITK_SOURCE_DIR}/Modules/Core/Common/include/itkFloatingPointExceptions.h")
if(NOT EXISTS "${_header}")
  message(FATAL_ERROR "Cannot find ITK header to patch: ${_header}")
endif()

file(READ "${_header}" _contents)
string(REPLACE "\n#include <cstdint>" "" _contents "${_contents}")
string(REPLACE
  "#include \"itkSingletonMacro.h\""
  "#include \"itkSingletonMacro.h\"\n#include <cstdint>"
  _contents
  "${_contents}")
file(WRITE "${_header}" "${_contents}")

set(_gdcm_mr3 "${ITK_SOURCE_DIR}/Modules/ThirdParty/GDCM/src/gdcm/Utilities/gdcmext/mec_mr3_io.c")
if(EXISTS "${_gdcm_mr3}")
  file(READ "${_gdcm_mr3}" _contents)
  if(NOT _contents MATCHES "#include <malloc.h>")
    string(REPLACE
      "#include <stdlib.h>"
      "#include <stdlib.h>\n#if defined(__MINGW32__)\n#include <malloc.h>\n#endif"
      _contents
      "${_contents}")
  endif()
  string(REPLACE
    "#ifdef _MSC_VER\n  return _aligned_malloc(size, alignment);"
    "#if defined(_MSC_VER) || defined(__MINGW32__)\n  return _aligned_malloc(size, alignment);"
    _contents
    "${_contents}")
  string(REPLACE
    "#ifdef _MSC_VER\n  _aligned_free(data.buffer);"
    "#if defined(_MSC_VER) || defined(__MINGW32__)\n  _aligned_free(data.buffer);"
    _contents
    "${_contents}")
  file(WRITE "${_gdcm_mr3}" "${_contents}")
endif()
