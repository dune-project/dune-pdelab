install(FILES assembler.hh
              coefficients.hh
              traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/operator/forward/instationary)

add_subdirectory(assembler)
