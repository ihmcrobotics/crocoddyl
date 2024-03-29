SET(${PROJECT_NAME}_EXAMPLES_PYTHON
  double_pendulum_fwddyn
  double_pendulum_invdyn
  quadrotor_fwddyn
  quadrotor_invdyn
  quadrotor_ubound
  arm_manipulation_fwddyn
  arm_manipulation_invdyn
  quadrupedal_gaits_fwddyn
  quadrupedal_gaits_invdyn
  quadrupedal_walk_ubound
  bipedal_walk_fwddyn
  bipedal_walk_invdyn
  bipedal_walk_ubound
  humanoid_manipulation
  humanoid_manipulation_ubound
  humanoid_taichi
  boxfddp_vs_boxddp
  )

INSTALL(FILES __init__.py DESTINATION ${PYTHON_SITELIB}/${PROJECT_NAME}/examples)
FOREACH(EXAMPLES ${${PROJECT_NAME}_EXAMPLES_PYTHON})
  PYTHON_BUILD(. "${EXAMPLES}.py")
  INSTALL(FILES ${EXAMPLES}.py DESTINATION ${PYTHON_SITELIB}/${PROJECT_NAME}/examples)
  ADD_CUSTOM_TARGET("examples-${EXAMPLES}"
    ${CMAKE_COMMAND} -E env PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} -c "import ${EXAMPLES}" \${INPUT})

  # examples are too slow in Debug mode to be used as tests
  IF(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    ADD_PYTHON_UNIT_TEST("example-python-${EXAMPLES}" "examples/${EXAMPLES}.py" bindings/python)
  ENDIF(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
ENDFOREACH(EXAMPLES ${${PROJECT_NAME}_EXAMPLES_PYTHON})

ADD_CUSTOM_TARGET("update_logfiles"
  ${CMAKE_COMMAND} -E env PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python:$ENV{PYTHONPATH} 
  PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} /bin/bash ${CMAKE_CURRENT_SOURCE_DIR}/log/update_logfiles.sh)
IF(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  IF(UNIX)
  ADD_TEST("test_logfiles"
${CMAKE_COMMAND} -E env PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python:$ENV{PYTHONPATH}
    PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} /bin/bash ${CMAKE_CURRENT_SOURCE_DIR}/log/update_logfiles.sh
    && /bin/bash ${CMAKE_CURRENT_SOURCE_DIR}/log/check_logfiles.sh)

    ADD_CUSTOM_TARGET("test_logfiles"
      ${CMAKE_COMMAND} -E env PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python:$ENV{PYTHONPATH}
      PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} /bin/bash ${CMAKE_CURRENT_SOURCE_DIR}/log/check_logfiles.sh)
  ENDIF(UNIX)
ENDIF(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
