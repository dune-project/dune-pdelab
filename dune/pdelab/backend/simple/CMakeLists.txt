set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/simple)
set(common_HEADERS
  descriptors.hh
  matrix.hh
  sparse.hh
  vector.hh)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
