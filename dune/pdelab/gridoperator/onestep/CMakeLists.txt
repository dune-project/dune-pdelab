install(FILES enginebase.hh
              jacobianengine.hh
              jacobianresidualengine.hh
              localassembler.hh
              nonlinearjacobianapplyengine.hh
              patternengine.hh
              prestageengine.hh
              residualengine.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/onestep)
