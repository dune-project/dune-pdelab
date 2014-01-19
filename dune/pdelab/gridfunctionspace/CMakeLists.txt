set(gridfunctionspacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridfunctionspace)
set(gridfunctionspace_HEADERS
  compositegridfunctionspace.hh
  genericdatahandle.hh
  gridfunctionspace.hh
  gridfunctionspacebase.hh
  gridfunctionspaceutilities.hh
  interpolate.hh
  intersectionindexset.hh
  localfunctionspace.hh
  localfunctionspacetags.hh
  localvector.hh
  powercompositegridfunctionspacebase.hh
  powergridfunctionspace.hh
  subspace.hh
  subspacelocalfunctionspace.hh
  tags.hh)

install(FILES ${gridfunctionspace_HEADERS} DESTINATION ${gridfunctionspacedir})
