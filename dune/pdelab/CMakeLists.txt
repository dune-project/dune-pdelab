set(pdelabdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab)
set(SUBDIRS                                        
        adaptivity                              
        backend                                 
        benchmark                               
        boilerplate                             
        common                                  
        constraints                             
        experimental                            
        finiteelement                   
        finiteelementmap                
        function                                
        gridfunctionspace               
        gridoperator                    
        gridoperatorspace               
        instationary                            
        linearsolver                    
        localoperator                   
        multistep                               
        newton                                  
        stationary                              
        test)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
