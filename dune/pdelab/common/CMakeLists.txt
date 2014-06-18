set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common)
set(common_HEADERS
  benchmarkhelper.hh
  borderindexidcache.hh
  clock.hh
  crossproduct.hh
  dofindex.hh
  elementmapper.hh
  exceptions.hh
  function.hh
  functionutilities.hh
  functionwrappers.hh
  geometrywrapper.hh
  globaldofindex.hh
  hostname.hh
  jacobiantocurl.hh
  logtag.hh
  multiindex.hh
  partitioninfoprovider.hh
  polymorphicbufferwrapper.hh
  range.hh
  simpledofindex.hh
  topologyutility.hh
  typetraits.hh
  unordered_map.hh
  unordered_set.hh
  utility.hh
  vtkexport.hh)

dune_add_library(dunepdelab
  clock.cc
  hostname.cc
  logtag.cc
  ADD_LIBS ${DUNE_LIBS})

install(FILES ${common_HEADERS} DESTINATION ${commondir})
