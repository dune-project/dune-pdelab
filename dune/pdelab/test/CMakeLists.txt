add_subdirectory(fem EXCLUDE_FROM_ALL)

pdelab_add_test(NAME testanalytic)
pdelab_add_test(NAME testclock)
pdelab_add_test(NAME testcomplexnumbers-istl SOURCES testcomplexnumbers.cc COMPILE_DEFINITIONS SUPERLU_NTYPE=3 USE_ISTL) # needed for complex SuperLU
pdelab_add_test(NAME testcomplexnumbers-istlbackend SOURCES testcomplexnumbers.cc COMPILE_DEFINITIONS SUPERLU_NTYPE=3 USE_ISTL_BACKEND) # needed for complex SuperLU
pdelab_add_test(NAME testconstraints COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME testfunction)
pdelab_add_test(NAME testgridfunctionspace)
pdelab_add_test(NAME testheat-instationary-periodic)
pdelab_add_test(NAME testlaplacedirichletccfv COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME testelasticity)
pdelab_add_test(NAME testlocalfunctionspace)
pdelab_add_test(NAME testpk2dinterpolation COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" ALBERTA_GRIDDIM 2)
pdelab_add_test(NAME testpk COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" ALBERTA_GRIDDIM 2)
pdelab_add_test(NAME testpoisson COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" ALBERTA_GRIDDIM 2)
pdelab_add_test(NAME testchunkedblockordering)
pdelab_add_test(NAME testrt0 COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" ALBERTA_GRIDDIM 2)
pdelab_add_test(NAME testrt02dgridfunctionspace COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" ALBERTA_GRIDDIM 2)
pdelab_add_test(NAME testutilities)
pdelab_add_test(NAME testinterpolate)
pdelab_add_test(NAME testrtfem COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME testbdmfem COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME testvectoriterator)
pdelab_add_test(NAME testpermutedordering)
pdelab_add_test(NAME testsimplebackend COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME test-transport-ccfv)

if(EIGEN3_FOUND)
  pdelab_add_test(NAME testeigenbackend COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
endif(EIGEN3_FOUND)

if(dune-alugrid_FOUND)
  pdelab_add_test(NAME testordering)
endif(dune-alugrid_FOUND)

if(GMP_FOUND)
  pdelab_add_test(NAME testopbfem)
  pdelab_add_test(NAME test-dg-amg)
endif(GMP_FOUND)

# Include periodic tests. Must run in parallel.
foreach(degree IN ITEMS 1 2)
  foreach(dg IN ITEMS 0 1)
    pdelab_add_test(NAME testpoisson-periodic-2d-deg${degree}-dg${dg}-parallel SOURCES testpoisson-periodic-2d.cc MPIRANKS 2 COMPILE_DEFINITIONS "DEGREE=${degree}" "DG=${dg}")
    pdelab_add_test(NAME testpoisson-periodic-2d-deg${degree}-dg${dg}-sequential SOURCES testpoisson-periodic-2d.cc COMPILE_DEFINITIONS "DEGREE=${degree}" "DG=${dg}")
    pdelab_add_test(NAME testpoisson-periodic-3d-deg${degree}-dg${dg}-parallel SOURCES testpoisson-periodic-3d.cc MPIRANKS 4 COMPILE_DEFINITIONS "DEGREE=${degree}" "DG=${dg}")
    pdelab_add_test(NAME testpoisson-periodic-3d-deg${degree}-dg${dg}-sequential SOURCES testpoisson-periodic-3d.cc COMPILE_DEFINITIONS "DEGREE=${degree}" "DG=${dg}")
  endforeach()
endforeach()

pdelab_add_test(NAME testnonoverlappingsinglephaseflow-yasp SOURCES testnonoverlappingsinglephaseflow.cc)
pdelab_add_test(NAME testnonoverlappingsinglephaseflow-yasp-parallel SOURCES testnonoverlappingsinglephaseflow.cc MPIRANKS 2)
if (UG_FOUND)
  pdelab_add_test(NAME testnonoverlappingsinglephaseflow-ug SOURCES testnonoverlappingsinglephaseflow.cc COMPILE_DEFINITIONS "GRID_UG")
  pdelab_add_test(NAME testnonoverlappingsinglephaseflow-ug-parallel SOURCES testnonoverlappingsinglephaseflow.cc COMPILE_DEFINITIONS "GRID_UG" MPIRANKS 2)
endif (UG_FOUND)

if (UG_FOUND)
  pdelab_add_test(NAME testnonoverlappingsinglephaseflow-boilerplate SOURCES testnonoverlappingsinglephaseflow-boilerplate.cc COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"" MPIRANKS 2)
endif (UG_FOUND)

if (UG_FOUND)
  pdelab_add_test(NAME testadaptivity)
  pdelab_add_test(NAME testnonoverlapping MPIRANKS 2)
endif (UG_FOUND)

set_tests_properties(
  testpk2dinterpolation
  testrt02dgridfunctionspace
  PROPERTIES SKIP_RETURN_CODE 77
)

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
