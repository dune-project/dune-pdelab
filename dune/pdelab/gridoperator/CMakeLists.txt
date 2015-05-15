add_subdirectory(common)
add_subdirectory(default)
add_subdirectory(onestep)

install(FILES gridoperator.hh
              onestep.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator)
