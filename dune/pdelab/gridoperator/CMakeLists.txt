add_subdirectory(common)
add_subdirectory(default)
add_subdirectory(onestep)
add_subdirectory(blockstructured)

install(FILES gridoperator.hh
              onestep.hh
              blockstructured.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator)
