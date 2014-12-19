add_subdirectory(fem EXCLUDE_FROM_ALL)

# TODO: install subfolder grids and all its contents
# TODO: install make_pvd.sh

# defined empty so we can add to it later
set(NORMALTESTS)
set(MOSTLYCLEANFILES)

set(noinst_HEADERS
        fmt.hh
        gnuplotgraph.hh
        gridexamples.hh
        l2difference.hh
        l2norm.hh)

set(GRIDDIM  2)

list(APPEND NORMALTESTS testanalytic)
add_executable(testanalytic testanalytic.cc)
target_link_libraries(testanalytic dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testanalytic GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES channel.vtu)

list(APPEND NORMALTESTS testclock)
add_executable(testclock testclock.cc)
target_link_libraries(testclock dunepdelab ${DUNE_LIBS})

# disable IntersectionIndexSet test for now, we know it's broken
# if(ALUGRID_FOUND)
#   list(APPEND NORMALTESTS test_composed_iis_gfs)
#   add_executable(test_composed_iis_gfs test-composed-iis-gfs.cc)
#   add_dune_alugrid_flags(test_composed_iis_gfs)
# endif(ALUGRID_FOUND )

list(APPEND NORMALTESTS testconstraints)
add_executable(testconstraints testconstraints.cc)
target_link_libraries(testconstraints dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testconstraints GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES testconstraintsp1.vtu testconstraintspowerp1.vtu)

list(APPEND NORMALTESTS testfunction)
add_executable(testfunction testfunction.cc)
target_link_libraries(testfunction dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testfunction GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES multi.vtu single.vtu)

list(APPEND NORMALTESTS testgridfunctionspace)
add_executable(testgridfunctionspace testgridfunctionspace.cc)
target_link_libraries(testgridfunctionspace dunepdelab ${DUNE_LIBS})

if(dune-alugrid_FOUND)
  list(APPEND NORMALTESTS testordering)
  add_executable(testordering testordering.cc)
  target_link_libraries(testordering dunepdelab ${DUNE_LIBS})

  list(APPEND NORMALTESTS testpermutedordering)
  add_executable(testpermutedordering testpermutedordering.cc)
  target_link_libraries(testpermutedordering dunepdelab ${DUNE_LIBS})
endif(ALUGRID_FOUND)
endif(dune-alugrid_FOUND)

list(APPEND NORMALTESTS testlaplacedirichletccfv)
add_executable(testlaplacedirichletccfv testlaplacedirichletccfv.cc)
target_link_libraries(testlaplacedirichletccfv dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testlaplacedirichletccfv)
add_dune_alberta_flags(testlaplacedirichletccfv GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES testlaplacedirichletccfv.vtu)

# This test still uses GridOperatorSpace
# list(APPEND NORMALTESTS testlaplacedirichletp12d)
# add_executable(testlaplacedirichletp12d testlaplacedirichletp12d.cc)
# target_link_libraries(testlaplacedirichletp12d dunepdelab ${DUNE_LIBS})
# add_dune_superlu_flags(testlaplacedirichletp12d)
# add_dune_alberta_flags(testlaplacedirichletp12d GRIDDIM 2)
# list(APPEND MOSTLYCLEANFILES testlaplacedirichletp12d.vtu)

list(APPEND NORMALTESTS testelasticity)
add_executable(testelasticity testelasticity.cc)
target_link_libraries(testelasticity dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testelasticity)
list(APPEND MOSTLYCLEANFILES testelasticity.vtu)

list(APPEND NORMALTESTS testlocalfunctionspace)
add_executable(testlocalfunctionspace testlocalfunctionspace.cc)
target_link_libraries(testlocalfunctionspace dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testlocalfunctionspace GRIDDIM 2)

list(APPEND NORMALTESTS testpk2dinterpolation)
add_executable(testpk2dinterpolation testpk2dinterpolation.cc)
target_link_libraries(testpk2dinterpolation dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testpk2dinterpolation GRIDDIM 2)

list(APPEND NORMALTESTS testpk)
add_executable(testpk testpk.cc)
target_link_libraries(testpk dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testpk GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES testpk.vtu)

list(APPEND NORMALTESTS testpoisson)
add_executable(testpoisson testpoisson.cc)
target_link_libraries(testpoisson dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testpoisson)
add_dune_alberta_flags(testpoisson GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES poisson_*.vtu)

list(APPEND NORMALTESTS testrt0)
add_executable(testrt0 testrt0.cc)
target_link_libraries(testrt0 dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testrt0 GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES testrt0.vtu)

list(APPEND NORMALTESTS testrt02dgridfunctionspace)
add_executable(testrt02dgridfunctionspace testrt02dgridfunctionspace.cc)
target_link_libraries(testrt02dgridfunctionspace dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testrt02dgridfunctionspace GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES rt02dgridfunctionspace-*.vtu)

list(APPEND NORMALTESTS testutilities)
add_executable(testutilities testutilities.cc)
target_link_libraries(testutilities dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testutilities GRIDDIM 2)
list(APPEND MOSTLYCLEANFILES edger.vtu interpolated.vtu q1.vtu taylorhood.vtu)

list(APPEND NORMALTESTS testinterpolate)
add_executable(testinterpolate testinterpolate.cc)
target_link_libraries(testinterpolate dunepdelab ${DUNE_LIBS})

if(PETSC_FOUND)
  list(APPEND NORMALTESTS testpoisson_petsc)
  add_executable(testpoisson_petsc testpoisson-petsc.cc)
  target_link_libraries(testpoisson_petsc dunepdelab ${DUNE_LIBS})
  add_dune_petsc_flags(testpoisson_petsc)
  list(APPEND MOSTLYCLEANFILES poisson_petsc_*.vtu)
endif(PETSC_FOUND )

if(GMP_FOUND)
  list(APPEND NORMALTESTS testopbfem)
  add_executable(testopbfem testopbfem.cc)
  target_link_libraries(testopbfem dunepdelab ${DUNE_LIBS})
  add_dune_gmp_flags(testopbfem)

  list(APPEND NORMALTESTS test-dg-amg)
  add_executable(test-dg-amg test-dg-amg.cc)
  target_link_libraries(test-dg-amg dunepdelab ${DUNE_LIBS})
  add_dune_gmp_flags(test-dg-amg)
endif(GMP_FOUND )

list(APPEND NORMALTESTS testrtfem)
add_executable(testrtfem testrtfem.cc)
target_link_libraries(testrtfem dunepdelab ${DUNE_LIBS})

list(APPEND NORMALTESTS testbdmfem)
add_executable(testbdmfem testbdmfem.cc)
target_link_libraries(testbdmfem dunepdelab ${DUNE_LIBS})

list(APPEND NORMALTESTS testvectoriterator)
add_executable(testvectoriterator testvectoriterator.cc)
target_link_libraries(testvectoriterator dunepdelab ${DUNE_LIBS})

list(APPEND NORMALTESTS testpermutedordering)
add_executable(testpermutedordering testpermutedordering.cc)
target_link_libraries(testpermutedordering dunepdelab ${DUNE_LIBS})

foreach(i ${NORMALTESTS})
  add_test(${i} ${i})
  # add different flags and directory to example grids
  set_property(TARGET ${i}
    APPEND PROPERTY
    COMPILE_DEFINITIONS "GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\"")
  add_dune_alugrid_flags(${i})
  add_dune_ug_flags(${i})
  add_dune_parmetis_flags(${i})
endforeach(i ${NORMALTESTS})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${NORMALTESTS})
