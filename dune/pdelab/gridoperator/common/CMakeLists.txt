exclude_from_headercheck(assembler.hh)

set(gridoperatorcommon_HEADERS
  assembler.hh
  assemblerutilities.hh
  borderdofexchanger.hh
  diagonallocalmatrix.hh
  gridoperatorutilities.hh
  localassemblerenginebase.hh
  localmatrix.hh
  timesteppingparameterinterface.hh
  )

install(FILES ${gridoperatorcommon_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/common)
