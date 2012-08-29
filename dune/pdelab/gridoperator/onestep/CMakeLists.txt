set(gridoperatoronestepdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/onestep)

set(gridoperatoronestep_HEADERS                    
        enginebase.hh                           
        jacobianengine.hh                       
        jacobianresidualengine.hh               
        localassembler.hh                       
        patternengine.hh                        
        prestageengine.hh                       
        residualengine.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperatoronestep_HEADERS} DESTINATION ${gridoperatoronestepdir})