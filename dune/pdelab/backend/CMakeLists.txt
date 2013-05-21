set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend)
set(common_HEADERS  backendselector.hh
                 eigenmatrixbackend.hh
                 eigensolverbackend.hh
                 eigenvectorbackend.hh
                 istlmatrixbackend.hh
                 istlsolverbackend.hh
                 istlvectorbackend.hh
                 novlpistlsolverbackend.hh
                 ovlpistlsolverbackend.hh
                 petscmatrixbackend.hh
                 petscnestedvectorbackend.hh
                 petscutility.hh
                 petscvectorbackend.hh
                 seqistlsolverbackend.hh
                 solver.hh
                 vectorutilities.hh)

add_subdirectory(istl)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})
