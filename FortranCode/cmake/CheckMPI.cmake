CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
FUNCTION(fortran_check_MPI)


find_program(MPIEXEC aprun)
find_package(MPI)

get_filename_component(MPIEXEC_NAME ${MPIEXEC} NAME)
get_filename_component(MPI_Fortran_COMPILER_NAME ${MPI_Fortran_COMPILER} NAME)

message(STATUS "MPI exec: ${MPIEXEC_NAME}")
message(STATUS "MPI compiler: ${MPI_Fortran_COMPILER_NAME}")

if (${MPIEXEC_NAME} STREQUAL "aprun" AND ${MPI_Fortran_COMPILER_NAME} STREQUAL "ftn")
  set(RUNNING_ON_CRAY TRUE)
else (${MPIEXEC_NAME} STREQUAL "aprun" AND ${MPI_Fortran_COMPILER_NAME} STREQUAL "ftn")
  set(RUNNING_ON_CRAY FALSE)
endif (${MPIEXEC_NAME} STREQUAL "aprun" AND ${MPI_Fortran_COMPILER_NAME} STREQUAL "ftn")

set(RUNNING_ON_CRAY ${RUNNING_ON_CRAY} CACHE STRING "We are running on a Cray system. Need aprun to launch Zacros")
message(STATUS "Are we running on Cray? ${RUNNING_ON_CRAY}")

ENDFUNCTION(fortran_check_MPI)
