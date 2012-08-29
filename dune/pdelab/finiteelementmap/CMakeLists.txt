set(finiteelementmapdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/finiteelementmap)

set(finiteelementmap_HEADERS                       
        conformingconstraints.hh                
        edges0.5fem.hh                          
        finiteelementmap.hh                     
        global.hh                               
        hangingnodeconstraints.hh               
        hangingnodemanager.hh                   
        l2orthonormal.hh                        
        mfdconstraints.hh                       
        mimeticfem.hh                           
        monomfem.hh                             
        opbfem.hh                               
        p0constraints.hh                        
        p0ghostconstraints.hh                   
        p0fem.hh                                
        p12dfem.hh                              
        p1fem.hh                                
        pk1dbasis.hh                            
        pk2dfem.hh                              
        pk3dfem.hh                              
        pkfem.hh                                
        powerfem.hh                             
        q12dfem.hh                              
        q1fem.hh                                
        q22dfem.hh                              
        qkdg.hh                                 
        rannacher_turek2dfem.hh                 
        rt02dfem.hh                             
        rt0constraints.hh                       
        rt0q2dfem.hh                            
        rt0q3dfem.hh                            
        rt0qfem.hh                              
        variablemonomfem.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${finiteelementmap_HEADERS} DESTINATION ${finiteelementmapdir})