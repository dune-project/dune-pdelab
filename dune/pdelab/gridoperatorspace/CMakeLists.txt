set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperatorspace)
set(common_HEADERS                                 
        gridoperatorspace.hh                    
        gridoperatorspaceutilities.hh           
        instationarygridoperatorspace.hh        
        localmatrix.hh                          
        printmatrix.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})