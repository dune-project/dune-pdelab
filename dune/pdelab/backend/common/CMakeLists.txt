set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/common)
set(common_HEADERS
  uncachedmatrixview.hh
  uncachedvectorview.hh)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
