if(NOT HAVE_EIGEN)
  exclude_from_headercheck(
    descriptors.hh
    matrix.hh
    solvers.hh
    vector.hh)
endif()

install(FILES descriptors.hh
              matrix.hh
              solvers.hh
              vector.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/eigen)
