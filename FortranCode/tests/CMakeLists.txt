if(dotests)
  # Regressions execute the program as a whole and then check the results
  function(add_regression NAME DATA LABEL MYTIMEOUT)
    if (CMAKE_CONFIGURATION_TYPES)
      add_test( NAME ${NAME}
                COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/sample_run.py
                        ${PROJECT_BINARY_DIR}/$<CONFIG>/${EXEC_NAME} ${DATA} ${ARGN} ${TEST_RUNNER_ARGS})
    else ()
      add_test( NAME ${NAME}
                COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/sample_run.py
                        ${PROJECT_BINARY_DIR}/${EXEC_NAME} ${DATA} ${ARGN} ${TEST_RUNNER_ARGS})
    endif ()
    set_tests_properties(${NAME} PROPERTIES LABELS "regression;${LABEL}" TIMEOUT ${MYTIMEOUT})
  endfunction(add_regression)

  # Fixtures are tests that start from a given set of input files, but then execute some routine of
  # their own
  function(add_fixture_test NAME DATA)
    file(RELATIVE_PATH test_type ${PROJECT_SOURCE_DIR}/tests/ ${CMAKE_CURRENT_SOURCE_DIR})
    set(exec_name test_${NAME})
    set(test_name ${test_type}_${NAME})
    add_executable(${exec_name} ${NAME}.f90)
    target_link_libraries(${exec_name} fixture)
    include_directories(${PROJECT_BINARY_DIR}/tests)
    if(doopenmp)
      set_target_properties(${exec_name} PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
    endif()

    if (CMAKE_CONFIGURATION_TYPES)
      add_test( NAME ${NAME}
                COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/sample_run.py
                        ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/${exec_name} ${DATA} ${TEST_RUNNER_ARGS}
                        --fixture
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    else ()
      add_test( NAME ${NAME}
                COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/sample_run.py
                        ${CMAKE_CURRENT_BINARY_DIR}/${exec_name} ${DATA} ${TEST_RUNNER_ARGS}
                        --fixture
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
	endif ()
    set_tests_properties(${NAME} PROPERTIES LABELS "${test_type};unittest;fast")
  endfunction(add_fixture_test)


  # Unit-tests are stand-alone programs that return non-zero on failure
  function(add_unittest NAME)
    file(RELATIVE_PATH test_type ${PROJECT_SOURCE_DIR}/tests/ ${CMAKE_CURRENT_SOURCE_DIR})
    set(exec_name test_${NAME})
    set(test_name ${test_type}_${NAME})
    add_executable(${exec_name} ${NAME}.f90 ${ARGN})
    target_link_libraries(${exec_name} meso_approx_lib)
    include_directories(${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/tests)
    if(doopenmp)
      set_target_properties(${exec_name} PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
    endif()

    add_test( NAME ${test_name} COMMAND ${exec_name}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    set_tests_properties(${test_name} PROPERTIES LABELS "${test_type};unittest;fast")
  endfunction(add_unittest)

  # Library to which fixture tests link
  include_directories(${PROJECT_BINARY_DIR})
  # add_library(fixture STATIC fixture.f90)
  # add_library(unit_tests STATIC unit_tests.f90)
  target_link_libraries(meso_approx_lib)

  set(fasttimeout   "100")
  set(mediumtimeout "1000")
  set(slowtimeout   "10000")

  add_regression( SAMPLE_RUN_BPEC                                   0 fast     ${fasttimeout})
  add_regression( SAMPLE_RUN_K2NNC2                                 1 fast     ${fasttimeout})
  add_regression( SAMPLE_RUN_K3NNC2                                 2 medium   ${mediumtimeout})

endif(dotests)

#add_subdirectory(caching)
#add_subdirectory(pattern_iterator)
