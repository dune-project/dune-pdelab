install(FILES adapter.hh
              line_search.hh
              operator.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/operator)

add_subdirectory(forward)
add_subdirectory(inverse)
add_subdirectory(local_assembly)
