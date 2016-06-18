add_subdirectory(common)
add_subdirectory(eigen)
add_subdirectory(istl)
add_subdirectory(simple)

if(NOT HAVE_EIGEN)
  exclude_from_headercheck(
    eigen.hh)
endif()

install(FILES eigen.hh
              interface.hh
              istl.hh
              simple.hh
              solver.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend)
