#install headers
install(FILES
  arpackpp_geneo.hh
  coarsespace.hh
  conbase_exterior.hh
  geneo.hh
  geneobasis.hh
  liptonbabuskabasis.hh
  localoperator_ovlp_region.hh
  multicommdatahandle.hh
  neumann_boundary_condition.hh
  onedimsubdomainbasis.hh
  partitionofunity.hh
  subdomainbasis.hh
  subdomainprojectedcoarsespace.hh
  two_level_schwarz.hh
  zembasis.hh

  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/backend/istl/geneo)
