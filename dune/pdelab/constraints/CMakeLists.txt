add_subdirectory(common)

install(FILES conforming.hh
              hangingnode.hh
              hangingnodemanager.hh
              interiornode.hh
              noconstraints.hh
              p0.hh
              p0ghost.hh
              raviartthomas0.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/constraints)
