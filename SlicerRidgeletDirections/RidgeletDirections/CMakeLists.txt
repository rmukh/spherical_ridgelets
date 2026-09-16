set(MODULE_NAME RidgeletDirections)

set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

set(MODULE_PYTHON_RESOURCES
  RidgeletDirectionsLib/__init__.py
  RidgeletDirectionsLib/core.py
  )

slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

# CLI Module for computing ridgelet directions
set(CLI_XML_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CLI)

slicerMacroBuildCLI(
  NAME RidgeletDirectionsCompute
  CLIENT_SOURCE_DIR ${CLI_XML_DIR}
  CLIENT_EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/CLI/RidgeletDirectionsCompute.py
  WITH_GENERIC_TESTS
  )