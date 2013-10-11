set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperatorspace)
set(common_HEADERS
    gridoperatorspaceutilities.hh
    localmatrix.hh
    printmatrix.hh)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
