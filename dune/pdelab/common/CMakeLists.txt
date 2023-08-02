if(NOT MINGW)
  set(clock_hh "clock.hh")
endif()

install(FILES algebra.hh
              benchmarkhelper.hh
              borderindexidcache.hh
              ${clock_hh}
              container_entry.hh
              container_traversal.hh
              crossproduct.hh
              dofindex.hh
              elementmapper.hh
              entity_cast.hh
              error_condition.hh
              exceptions.hh
              for_each.hh
              function.hh
              functionutilities.hh
              functionwrappers.hh
              geometrywrapper.hh
              globaldofindex.hh
              hostname.hh
              instationaryfilenamehelper.hh
              intersectiontype.hh
              jacobiantocurl.hh
              local_container_entry.hh
              local_container.hh
              local_matrix.hh
              logtag.hh
              multiindex.hh
              partitionviewentityset.hh
              polymorphicbufferwrapper.hh
              property_tree.hh
              property.hh
              quadraturerules.hh
              range.hh
              referenceelements.hh
              simpledofindex.hh
              topologyutility.hh
              trace.hh
              tree_traversal.hh
              typetraits.hh
              utility.hh
              vtkexport.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common)

add_subdirectory(concurrency)
add_subdirectory(convergence)
