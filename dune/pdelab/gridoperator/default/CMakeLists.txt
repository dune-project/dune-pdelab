set(gridoperatordefaultdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/default)

set(gridoperatordefault_HEADERS                            
        assembler.hh                                    
        jacobianengine.hh
        jacobianapplyengine.hh
        localassembler.hh                               
        patternengine.hh                                
        residualengine.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperatordefault_HEADERS} DESTINATION ${gridoperatordefaultdir})