set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/dense)
set(common_HEADERS
  descriptors.hh
  matrix.hh
  vector.hh)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
