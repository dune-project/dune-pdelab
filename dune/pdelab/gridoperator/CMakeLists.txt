set(gridoperatordir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator)
set(SUBDIRS                                        
        common                                  
        default                                 
        onestep
        tbb)

set(gridoperator_HEADERS                           
        gridoperator.hh                         
        onestep.hh
        tbb.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperator_HEADERS} DESTINATION ${gridoperatordir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
