add_subdirectory(common)
add_subdirectory(default)
add_subdirectory(fastdg)
add_subdirectory(onestep)

install(FILES fastdg.hh
              gridoperator.hh
              onestep.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator)
