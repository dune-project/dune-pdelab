# This module defines a macro to add tests in pdelab.
# Adding a new test is a one-liner.
#
# pdelab_add_test(NAME testname
#                [SOURCES src1 [, src2, ...]]
#                [COMMAND cmd [args]]
#                [MPIRANKS ranks]
#                [COMPILE_DEFINITIONS def1 [, def2, ...]]
#                [COMPILE_OPTIONS opt1 [, opt2, ...]]
#                [ALBERTA_GRIDDIM gdim]
#                [ALBERTA_WORLDDIM wdim]
#  )
#
# The macro will do the following steps:
#  * add an executable called testname, that depends on the given source files. If no source files are given, the file <testname>.cc is considered a dependency.
#  * sets additional flags and defines on the target (this allows to write loops over lists of flags)
#  * registers the test
#
# The following features can be used optionally:
#  * parallel test execution.
#  * add flags for the alberta grid manager through ALBERTA_{GRID,WORLD}DIM. This is necessary,
#    as Alberta is the only external package that cannot be handled through dune_enable_all_packages()

# This target will be used to build all tests
add_custom_target(build_tests)

function(pdelab_add_test)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS NAME MPIRANKS ALBERTA_GRIDDIM ALBERTA_WORLDDIM)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_OPTIONS COMMAND)
  cmake_parse_arguments(PDELABTEST "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if(PDELABTEST_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${PDELABTEST_UNPARSED_ARGUMENTS}') for pdelab_add_test()!")
  endif()

  # by default, a test is built from a file with the same name and a ".cc" suffix
  if("${PDELABTEST_SOURCES}" STREQUAL "")
    set(PDELABTEST_SOURCES "${PDELABTEST_NAME}.cc")
  endif()

  add_executable(
    "${PDELABTEST_NAME}"
    EXCLUDE_FROM_ALL
    ${PDELABTEST_SOURCES}
    )

  add_dependencies(build_tests "${PDELABTEST_NAME}")

  target_compile_definitions(
    "${PDELABTEST_NAME}"
    PUBLIC
    ${PDELABTEST_COMPILE_DEFINITIONS}
    )

  target_compile_options(
    "${PDELABTEST_NAME}"
    PUBLIC
    ${PDELABTEST_COMPILE_OPTIONS}
    )

  if(NOT "${PDELABTEST_ALBERTA_GRIDDIM}" STREQUAL "")
    if("${PDELABTEST_ALBERTA_WORLDDIM}" STREQUAL "")
      set(PDELABTEST_ALBERTA_WORLDDIM ${PDELABTEST_ALBERTA_GRIDDIM})
    endif()
    add_dune_alberta_flags(${PDELABTEST_NAME} GRIDDIM ${PDELABTEST_ALBERTA_GRIDDIM} WORLDDIM ${PDELABTEST_ALBERTA_WORLDDIM})
  endif()

  if("${PDELABTEST_COMMAND}" STREQUAL "")
    set(PDELABTEST_COMMAND "${PDELABTEST_NAME}")
  endif()

  set(register_test TRUE)

  if(NOT "${PDELABTEST_MPIRANKS}" STREQUAL "")
    if(NOT "${PDELABTEST_MPIRANKS}" MATCHES "[1-9][0-9]*")
      message(ERROR "While creating test '${PDELABTEST_NAME}: invalid number of MPI ranks (${PDELABTEST_MPIRANKS})")
    endif()

    if(MPI_FOUND)
      set(PDELABTEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PDELABTEST_MPIRANKS} ${MPIEXEC_PREFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/${PDELABTEST_COMMAND}" ${MPIEXEC_POSTFLAGS})
    else()
      message(WARNING "Test '${PDELABTEST_NAME}' requires MPI, but MPI was not found. Test will be built, but not run")
      set(register_test FALSE)
    endif()
  endif()

  if(${register_test})
    # by default, the test is run by simply invoking the built executable
    _add_test(
      NAME ${PDELABTEST_NAME}
      COMMAND ${PDELABTEST_COMMAND}
      )
  endif()
endfunction()

# Override the builtin add_test command to give a warning if used from within dune-pdelab
macro(add_test)
  if(CMAKE_PROJECT_NAME STREQUAL dune-pdelab)
    message(WARNING "You are using the command add_test from within dune-pdelab. Please use pdelab_add_test instead.")
  endif()
  _add_test(${ARGN})
endmacro()
