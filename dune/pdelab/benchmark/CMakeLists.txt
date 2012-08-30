# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
set(NORMALTESTS  )
set(MOSTLYCLEANFILES )
set(check_SCRIPTS  )

# programs just to build when "make check" is used
set(check_PROGRAMS  ${NORMALTESTS})

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# TRUNK TARGETS FOLLOW

list(APPEND NORMALTESTS poisson-old)
add_executable(poisson-old poisson-old.cc)
list(APPEND MOSTLYCLEANFILES poisson-old_*.vtu poisson-old_*timings.txt)

list(APPEND NORMALTESTS laplacedirichletccfv-old)
add_executable(laplacedirichletccfv-old laplacedirichletccfv-old.cc)
list(APPEND MOSTLYCLEANFILES laplacedirichletccfv-old_*.vtu laplacedirichletccfv-old_*timings.txt)

foreach(_test ${NORMALTESTS})
  add_dune_ug_flags(${_test})
  add_dune_alugrid_flags(${_test})
  add_dune_amiramesh_flags(${_test})
  add_dune_mpi_flags(${_test})
  add_dune_parmetis_flags(${_test})
endforeach(_test ${NORMALTEST})

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

set(DISTCHECK_CONFIGURE_FLAGS  --with-dune=${DUNEROOT} CXX="${CXX}" CC="${CC}")

# list of tests to run
set(TESTS  ${NORMALTESTS})

add_test(NAME poisson-old WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/poisson-old ${CMAKE_CURRENT_SOURCE_DIR}/poisson-old.ini)
add_test(NAME laplacedirichletccfv-old WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/laplacedirichletccfv-old ${CMAKE_CURRENT_SOURCE_DIR}/laplacedirichletccfv-old.ini)


# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})