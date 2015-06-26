add_subdirectory(fem EXCLUDE_FROM_ALL)

pdelab_add_test(NAME testanalytic)
pdelab_add_test(NAME testclock)
pdelab_add_test(NAME testconstraints COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
pdelab_add_test(NAME testfunction)
pdelab_add_test(NAME testgridfunctionspace)
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

if(dune-alugrid_FOUND)
  pdelab_add_test(NAME testordering)
endif(dune-alugrid_FOUND)

if(GMP_FOUND)
  pdelab_add_test(NAME testopbfem)
  pdelab_add_test(NAME test-dg-amg)
endif(GMP_FOUND)

set_tests_properties(
  testpk2dinterpolation
  testrt02dgridfunctionspace
  PROPERTIES SKIP_RETURN_CODE 77
)
