install(FILES
      identity.hh
      iteratorsplit.hh
      metis.hh
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common/partition)

add_subdirectory(halo)
