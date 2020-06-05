add_subdirectory(fem)

dune_add_test(SOURCES testcombinedoperator.cc)

dune_add_test(SOURCES testconvectiondiffusiondg.cc
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(SOURCES testnewton.cc
              CMAKE_GUARD SUPERLU_FOUND)

dune_add_test(SOURCES testoldnewton.cc)

dune_add_test(SOURCES testinstationary.cc)

dune_add_test(SOURCES testanalytic.cc)

dune_add_test(SOURCES testbindtime.cc)

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
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpk2dinterpolation)

dune_add_test(SOURCES testpk.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpk)

dune_add_test(SOURCES testpoisson.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testpoisson)

dune_add_test(SOURCES testmatrixfree.cc)

dune_add_test(SOURCES testchunkedblockordering.cc)

dune_add_test(SOURCES testrt0.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testrt0)

dune_add_test(SOURCES testrt02dgridfunctionspace.cc
              COMPILE_DEFINITIONS GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")
add_dune_alberta_flags(testrt02dgridfunctionspace)

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

dune_add_test(SOURCES testordering.cc
              CMAKE_GUARD dune-alugrid_FOUND)

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

exclude_from_headercheck(
  gridexamples.hh
  helmholtzoperator.hh
  testadaptivity-adaptivity.hh
  testadaptivity-bcextension.hh
  testadaptivity-bctype.hh
  testadaptivity-error-indicator.hh
  testadaptivity-operator.hh
  testcomplexnumbers-problem.hh
  testnonoverlappingsinglephaseflow-boilerplate-problem.hh
  testnonoverlappingsinglephaseflow-problem.hh
  testpoisson-periodic.hh
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vtk)
