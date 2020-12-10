CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project ("MesoscopicKineticApproximations" Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
option(doopenmp    "Will compile with openmp, if found." off)

include(compilerinfo2008_Fortran)
include(CheckMPI)
include(VersionAndGitRef)

set_version(0.01dev)
get_gitref()

set(F90_SOURCES fdjac.f90
				fminln.f90
				func.f90
				funcv.f90
				global_constants.f90
				init_1NN.f90
				init_2NN.f90
				init_3NN.f90
				lnsrch.f90
				lubksb.f90
				ludcmp.f90
				meso_approx.f90
				meso_approx_inst.f90
				newt.f90
				nr.f90
				nrtype.f90
				nrutil.f90
				parser_module.f90
				constants_module.f90
				err_warn_module.F90
				calculation_setup_module.f90
				calculation_parser_module.f90
				polint.f90
				powell.f90)

set(FREELINE_SOURCES fdjac.f90
				fminln.f90
				func.f90
				funcv.f90
				global_constants.f90
				init_1NN.f90
				init_2NN.f90
				init_3NN.f90
				lnsrch.f90
				lubksb.f90
				ludcmp.f90
				meso_approx.f90
				meso_approx_inst.f90
				newt.f90
				nr.f90
				nrtype.f90
				nrutil.f90
				parser_module.f90
				polint.f90
				powell.f90)

enable_language(Fortran)

fortran_check_compiler2008()

# configure_file(info_module.F90 info_module_processed.F90)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  SET_SOURCE_FILES_PROPERTIES( ${FREELINE_SOURCES} PROPERTIES COMPILE_FLAGS
    -ffree-line-length-0 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -kind=byte -ieee=full")
endif()


if(doopenmp)
  find_package(OpenMP_Fortran)
  if(OpenMP_Fortran_FLAGS)
    set(OPENMP_SOURCES  meso_approx.f90 main_meso_approx.f90) # placeholder - not used currently
    foreach(OPENMP_SOURCE ${OPENMP_SOURCES})
      get_source_file_property(OLD_FLAGS ${OPENMP_SOURCE} COMPILE_FLAGS)
      if (OLD_FLAGS)
        MESSAGE(STATUS "The old flags are ${OLD_FLAGS} on ${OPENMP_SOURCE}")
        SET(All_FLAGS "${OLD_FLAGS} ${OpenMP_Fortran_FLAGS}")
      else (OLD_FLAGS)
        SET(All_FLAGS "${OpenMP_Fortran_FLAGS}")
      endif (OLD_FLAGS)
      set_source_files_properties(${OPENMP_SOURCE}
                                   PROPERTIES
                                   COMPILE_FLAGS ${All_FLAGS})
    endforeach(OPENMP_SOURCE)
  else(OpenMP_Fortran_FLAGS)
    message("No OpenMP flag found.")
    set(doopenmp FALSE)
  endif(OpenMP_Fortran_FLAGS)
endif(doopenmp)

fortran_check_MPI()

if (MSVC OR MSYS OR MINGW)
  message(STATUS "Compiling for Windows platform (MS Visual Studio, MSys, MinGW, ...)")
  set (EXEC_NAME "MesoApprox")
elseif (WIN32)
  message(STATUS "Compiling for Windows platform")
  set (EXEC_NAME "MesoApprox")
elseif (UNIX AND NOT APPLE)
  message(STATUS "Compiling for Linux platform (BSD, Solaris, Minix, ...)")
  set (EXEC_NAME "mesoapprox.x")
elseif (APPLE)
  message(STATUS "Compiling for MacOS platform")
  set (EXEC_NAME "mesoapprox.x")
else ()
  set (EXEC_NAME "mesoapprox.x")
endif ()

if(CMAKE_CONFIGURATION_TYPES)
    message("Multi-configuration generator")
else()
    message("Single-configuration generator")
endif()

add_library(meso_approx_lib STATIC ${F90_SOURCES})
add_executable(${EXEC_NAME} main_meso_approx.f90)
target_link_libraries(${EXEC_NAME} meso_approx_lib)
if(doopenmp)
  set_target_properties(${EXEC_NAME} PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()

find_package(PythonInterp)

if (RUNNING_ON_CRAY)
  SET(TEST_RUNNER_ARGS "--execrunner;aprun")
else (RUNNING_ON_CRAY)
  SET(TEST_RUNNER_ARGS )
endif (RUNNING_ON_CRAY)

SET(TEST_RUNNER_ARGS ${TEST_RUNNER_ARGS} CACHE STRING "MPI type runner passed to sample_run.py in tests")
message(STATUS "Test runner is ${TEST_RUNNER_ARGS}")

if(PYTHONINTERP_FOUND AND PYTHON_VERSION_STRING VERSION_GREATER "2.7.0" )
  set(dotests TRUE)
  message(STATUS "Enabling tests")
  # Use include test to suppress warning about DART
  # We are not actually using DART but turn it on to get the XML files that
  # can be converted to JUNIT for Jenkins
  include(CTest)
  ENABLE_TESTING()
  add_subdirectory(tests)
else(PYTHONINTERP_FOUND AND PYTHON_VERSION_STRING VERSION_GREATER "2.7.0")
  if(PYTHONINTERP_FOUND)
    message(STATUS "At least Python 2.7 is needed, Found ${PYTHON_VERSION_STRING}. Cannot run tests.")
  else(PYTHONINTERP_FOUND)
    message(STATUS "Python not found. Cannot run tests.")
  endif(PYTHONINTERP_FOUND)
  
  set(dotests FALSE)
endif(PYTHONINTERP_FOUND AND PYTHON_VERSION_STRING VERSION_GREATER "2.7.0" )

install (TARGETS ${EXEC_NAME} DESTINATION bin)

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#   message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
