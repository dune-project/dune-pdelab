set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/eigen)
set(common_HEADERS
  descriptors.hh
  matrix.hh
  solvers.hh
  vector.hh)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
