#install headers
install(FILES
  arpackpp_geneo.hh
  coarsespace.hh
  geneo.hh
  geneobasis.hh
  liptonbabuskabasis.hh
  localoperator_ovlp_region.hh
  multicommdatahandle.hh
  partitionofunity.hh
  subdomainbasis.hh
  subdomainprojectedcoarsespace.hh
  two_level_schwarz.hh

  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/istl/geneo)

install(FILES
  communicator_with_rank.hh
  geneobasis.hh
  geneo_matrix_setup.hh
  novlp_geneo_preconditioner.hh
  novlp_operators.hh
  overlaptools.hh
  schwarznonoverlapping.hh
  subdomainprojectedcoarsespace.hh
  two_level_schwarz.hh
  variablesizecommunicator_with_rank.hh

  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/istl/geneo/nonoverlapping)
