install(FILES basis.hh
              discrete_global_function.hh
              merging_strategy.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/basis)

add_subdirectory(backend)
add_subdirectory(constraints)
add_subdirectory(ordering)
add_subdirectory(prebasis)
add_subdirectory(wrapper)
