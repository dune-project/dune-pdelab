set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/eigen)
set(common_HEADERS
  descriptors.hh
  matrix.hh
  solvers.hh
  vector.hh)

if(NOT HAVE_EIGEN)
  exclude_from_headercheck(
    descriptors.hh
    matrix.hh
    solvers.hh
    vector.hh)
endif()

install(FILES ${common_HEADERS} DESTINATION ${commondir})
