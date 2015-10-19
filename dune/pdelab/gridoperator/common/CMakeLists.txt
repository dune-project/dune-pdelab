exclude_from_headercheck(assembler.hh)

install(FILES assembler.hh
              assemblerutilities.hh
              borderdofexchanger.hh
              diagonallocalmatrix.hh
              gridoperatorutilities.hh
              localassemblerenginebase.hh
              localmatrix.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/common)
