set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/multistep)
set(my_HEADERS                                     
        cache.hh                                
        gridoperatorspace.hh                    
        method.hh                               
        parameter.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${my_HEADERS} DESTINATION ${mydir})