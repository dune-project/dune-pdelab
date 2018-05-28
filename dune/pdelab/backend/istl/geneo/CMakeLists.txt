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
