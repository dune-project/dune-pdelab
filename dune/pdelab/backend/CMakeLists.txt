add_subdirectory(common)
add_subdirectory(eigen)
add_subdirectory(istl)
add_subdirectory(simple)

if(NOT HAVE_EIGEN)
  exclude_from_headercheck(
    eigen.hh
    eigenmatrixbackend.hh
    eigensolverbackend.hh
    eigenvectorbackend.hh)
endif()

install(FILES backendselector.hh
              eigen.hh
              eigenmatrixbackend.hh
              eigensolverbackend.hh
              eigenvectorbackend.hh
              istlmatrixbackend.hh
              istlsolverbackend.hh
              istlvectorbackend.hh
              novlpistlsolverbackend.hh
              ovlpistlsolverbackend.hh
              seqistlsolverbackend.hh
              simple.hh
              solver.hh
              tags.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend)
