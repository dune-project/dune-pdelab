add_subdirectory(common)
add_subdirectory(default)
add_subdirectory(onestep)
add_subdirectory(tbb)

install(FILES gridoperator.hh
              onestep.hh
              tbb.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator)
