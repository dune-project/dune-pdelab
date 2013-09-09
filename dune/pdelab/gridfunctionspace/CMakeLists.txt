set(gridfunctionspacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridfunctionspace)
set(gridfunctionspace_HEADERS
  compositegridfunctionspace.hh
  compositeorderingutilities.hh
  debug.hh
  dofinfo.hh
  genericdatahandle.hh
  gridfunctionspace.hh
  gridfunctionspacebase.hh
  gridfunctionspaceutilities.hh
  interpolate.hh
  intersectionindexset.hh
  leafordering.hh
  lexicographicordering.hh
  localfunctionspace.hh
  localfunctionspacetags.hh
  localvector.hh
  nonleaforderingbase.hh
  orderingbase.hh
  orderinginterface.hh
  powercompositegridfunctionspacebase.hh
  powergridfunctionspace.hh
  printvector.hh
  subspace.hh
  subspacelocalfunctionspace.hh
  tags.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridfunctionspace_HEADERS} DESTINATION ${gridfunctionspacedir})