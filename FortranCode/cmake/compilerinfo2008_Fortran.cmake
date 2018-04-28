CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

FUNCTION(fortran_check_compiler2008)
    message(STATUS "Checking if compiler supports compiler_version, compiler_options from Fortran 2008.")
    FILE (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortracompiler2008.f90"
    "
    program testcompiler2008
        use iso_fortran_env, only: compiler_version, compiler_options
        implicit none
        write(*,*) compiler_version()
        write(*,*) compiler_options()
    end program testcompiler2008
    ")

    TRY_RUN(COMPILER2008_RUN_RESULT COMPILER2008_COMPILE_RESULT ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortracompiler2008.f90
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        CMAKE_FLAGS ${CMAKE_REQUIRED_FLAGS}
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
        RUN_OUTPUT_VARIABLE RUN_OUTPUT)
    add_definitions(-DUSING_CMAKE)
    # Helper to know if we are using cmake as an alternative to fortran 2008
    IF (COMPILER2008_COMPILE_RESULT AND NOT COMPILER2008_RUN_RESULT)
        # COMPILER2008_RUN_RESULT is an exit code = 0 if the code runs
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                 "Determining if the Fortran compiler supports compiler_version and compiler_info passed with "
                 "the following output:\n${COMPILE_OUTPUT}\n\n")
        message(STATUS "Compiler supports compiler_versions and compiler_options.")
        SET (COMPILER_SUPPORTS_COMPILER2008_INTERNAL 1)
        add_definitions(-DCOMPILER_SUPPORTS_COMPILER2008)
    ELSE (COMPILER2008_COMPILE_RESULT AND NOT COMPILER2008_RUN_RESULT)
        FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                 "Determining if the Fortran compiler supports compiler_version and compiler_info failed with "
                 "the following output:\n${COMPILE_OUTPUT}\n\n")
        message(STATUS "Compiler does not support compiler_version and compiler_info.")
        SET (COMPILER_SUPPORTS_COMPILER2008_INTERNAL 0)
    ENDIF (COMPILER2008_COMPILE_RESULT AND NOT COMPILER2008_RUN_RESULT)

    SET (USING_CMAKE 1 CACHE STRING "Compiling with CMake")
    SET (COMPILER_SUPPORTS_COMPILER2008 ${COMPILER_SUPPORTS_COMPILER2008_INTERNAL}
         CACHE STRING "Fortran compiler supports compiler_version and compiler_info")
    MARK_AS_ADVANCED(COMPILER_SUPPORTS_COMPILER2008)


ENDFUNCTION(fortran_check_compiler2008)
