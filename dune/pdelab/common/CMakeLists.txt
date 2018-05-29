if(NOT MINGW)
  set(clock_hh "clock.hh")
endif()

install(FILES benchmarkhelper.hh
              borderindexidcache.hh
              ${clock_hh}
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
              instationaryfilenamehelper.hh
              intersectiontype.hh
              jacobiantocurl.hh
              logtag.hh
              multiindex.hh
              partitionviewentityset.hh
              polymorphicbufferwrapper.hh
              quadraturerules.hh
              range.hh
              referenceelements.hh
              simpledofindex.hh
              topologyutility.hh
              typetraits.hh
              utility.hh
              vtkexport.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common)
