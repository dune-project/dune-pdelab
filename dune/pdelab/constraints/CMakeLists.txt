set(constraintsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/constraints)
set(constraints_HEADERS
    conforming.hh
    hangingnode.hh
    hangingnodemanager.hh
    interiornode.hh
    interiornodeconstraints.hh
    mfd.hh
    p0.hh
    p0ghost.hh
    raviartthomas0.hh)

add_subdirectory(common)

install(FILES ${constraints_HEADERS} DESTINATION ${constraintsdir})
