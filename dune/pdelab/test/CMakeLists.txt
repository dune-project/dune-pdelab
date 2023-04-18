project(dune-pdelab-tests CXX)
cmake_minimum_required(VERSION 3.6)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-pdelab as if this was a downstream project
  find_package(dune-pdelab REQUIRED IMPORTED)
  dune_project()
endif()

# configuration options
option(COVERAGE_REPORT  "Enable coverage compiler flags")

# set up test library
enable_testing()

find_package(GTest)

if (NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.13.0
  )

  FetchContent_MakeAvailable(googletest)
  include(GoogleTest)

  set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES EXCLUDE_FROM_ALL TRUE)

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
    set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)
  endif()
endif()

# create main library target
dune_add_library(dunepdelab-test INTERFACE NO_EXPORT)

target_compile_options(dunepdelab-test INTERFACE
  -Wall
  -pedantic
  $<$<BOOL:${COVERAGE_REPORT}>:--coverage>
)

target_link_libraries(dunepdelab-test INTERFACE
  dunepdelab
  $<$<BOOL:${COVERAGE_REPORT}>:gcov>
)

# include tests executables
add_subdirectory(fem)
add_subdirectory(matrixfree)

dune_add_test(SOURCES testcombinedoperator.cc)

dune_add_test(SOURCES testchecklopinterface.cc)

dune_add_test(SOURCES testconvectiondiffusiondg.cc
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(SOURCES testnewton.cc
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(SOURCES testoldnewton.cc)

dune_add_test(SOURCES testinstationary.cc)

dune_add_test(SOURCES testanalytic.cc)

dune_add_test(SOURCES testbindtime.cc)

dune_add_test(SOURCES testblocklocaloperatorwrapper.cc)

dune_add_test(NAME testcomplexnumbers-istl
              SOURCES testcomplexnumbers.cc
              COMPILE_DEFINITIONS SUPERLU_NTYPE=3
                                  USE_ISTL
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(NAME testcomplexnumbers-istlbackend
              SOURCES testcomplexnumbers.cc
              COMPILE_DEFINITIONS SUPERLU_NTYPE=3
                                  USE_ISTL_BACKEND
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(SOURCES testconstraints.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")

dune_add_test(SOURCES testdunefunctionsgfs.cc
              MPI_RANKS 1 2 4
              TIMEOUT 300)

dune_add_test(SOURCES testelectrodynamic.cc
              CMD_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/grids/circ_in_rect.msh"
              CMAKE_GUARD dune-alugrid_FOUND SUPERLU_FOUND)

dune_add_test(SOURCES testfunction.cc)

dune_add_test(SOURCES testgridfunctionspace.cc)

dune_add_test(SOURCES testgridfunctionspace-fixedsize.cc)

dune_add_test(SOURCES testpowergridfunctionspace-fixedsize.cc)

dune_add_test(SOURCES testplasticitygfs.cc)

dune_add_test(SOURCES testheat-instationary-periodic.cc)

dune_add_test(SOURCES testelasticity.cc)

dune_add_test(SOURCES testgeneo.cc
              MPI_RANKS 2
              TIMEOUT 300
              CMAKE_GUARD SuiteSparse_UMFPACK_FOUND ARPACKPP_FOUND
              )

dune_add_test(SOURCES testdglegendre.cc)

dune_add_test(SOURCES testfastdgassembler.cc)

dune_add_test(SOURCES testinstationaryfastdgassembler.cc)

dune_add_test(SOURCES testlocalfunctionspace.cc)

dune_add_test(SOURCES testlocalmatrix.cc)

dune_add_test(SOURCES testlocaloperatorinterface.cc)

dune_add_test(SOURCES testpk2dinterpolation.cc
              TIMEOUT 600
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpk2dinterpolation WORLDDIM 2)

dune_add_test(SOURCES testpk.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpk WORLDDIM 2)

dune_add_test(SOURCES testpoisson.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpoisson WORLDDIM 2)

dune_add_test(SOURCES testpoisson-multidomain.cc
              COMPILE_DEFINITIONS
                GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"
                DUNE_PDELAB_ENABLE_EXPERIMENTAL_MULTIDOMAIN_SUPPORT=1
              CMAKE_GUARD dune-multidomaingrid_FOUND)

dune_add_test(SOURCES testmatrixfree.cc)

dune_add_test(SOURCES testchunkedblockordering.cc)

dune_add_test(SOURCES testrt0.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testrt0 WORLDDIM 2)

dune_add_test(SOURCES testrt02dgridfunctionspace.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testrt02dgridfunctionspace WORLDDIM 2)

dune_add_test(SOURCES testutilities.cc)

dune_add_test(SOURCES testinterpolate.cc)

dune_add_test(SOURCES testnumericaljacobianmethods.cc)

dune_add_test(SOURCES testrtfem.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")

dune_add_test(SOURCES testbdmfem.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")

dune_add_test(SOURCES testvectoriterator.cc)

dune_add_test(SOURCES testpermutedordering.cc)

dune_add_test(SOURCES testsimplebackend.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")

dune_add_test(SOURCES testclock.cc
              CMAKE_GUARD "NOT MINGW")

dune_add_test(SOURCES test-instationary-with-boundary-constraints.cc
              CMAKE_GUARD "NOT MINGW")

dune_add_test(SOURCES test-transport-ccfv.cc
              CMAKE_GUARD "NOT MINGW")

dune_add_test(SOURCES testeigenbackend.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"
              CMAKE_GUARD EIGEN3_FOUND)

dune_add_test(SOURCES testordering.cc)

dune_add_test(SOURCES testloadbalancing.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"
              MPI_RANKS 2
              TIMEOUT 300
              CMAKE_GUARD dune-alugrid_FOUND)

dune_add_test(SOURCES testopbfem.cc
              CMAKE_GUARD GMP_FOUND)

dune_add_test(SOURCES test-dg-amg.cc
              CMAKE_GUARD GMP_FOUND)

# Include periodic tests. Must run in parallel.
foreach(degree IN ITEMS 1 2)
  foreach(dg IN ITEMS 0 1)
    dune_add_test(NAME testpoisson-periodic-2d-deg${degree}-dg${dg}-parallel
                  SOURCES testpoisson-periodic-2d.cc
                  MPI_RANKS 1 2
                  COMPILE_DEFINITIONS DEGREE=${degree} DG=${dg}
                  TIMEOUT 300)

    dune_add_test(NAME testpoisson-periodic-3d-deg${degree}-dg${dg}-parallel
                  SOURCES testpoisson-periodic-3d.cc
                  MPI_RANKS 1 4
                  COMPILE_DEFINITIONS DEGREE=${degree} DG=${dg}
                  TIMEOUT 300)
  endforeach()
endforeach()

dune_add_test(NAME testnonoverlappingsinglephaseflow-yasp
              SOURCES testnonoverlappingsinglephaseflow.cc
              MPI_RANKS 1 2
              TIMEOUT 300)

dune_add_test(NAME testnonoverlappingsinglephaseflow-ug
              SOURCES testnonoverlappingsinglephaseflow.cc
              COMPILE_DEFINITIONS GRID_UG
              MPI_RANKS 1 2
              TIMEOUT 300
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES testnonoverlappingsinglephaseflow-boilerplate.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"
              MPI_RANKS 2
              TIMEOUT 300
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES testadaptivity.cc
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES testnonoverlapping.cc
              MPI_RANKS 2
              TIMEOUT 300
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES testranktransmission.cc
              MPI_RANKS 2
              TIMEOUT 300)

dune_add_test(SOURCES testl2.cc)

dune_add_test(SOURCES test-blocked-istl-ordering.cc)

dune_add_test(SOURCES testtimedependentboundary_ovlpqk.cc)

dune_add_test(SOURCES testnewton-printparameters.cc)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vtk)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-pdelab as if this was a downstream project
  finalize_dune_project()
endif()
