set(typetreedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common/typetree)

set(typetree_HEADERS                               
        accumulate_static.hh                    
        applytochildrensingletree.hh            
        applytochildrentreepair.hh              
        childextraction.hh                      
        compositenode.hh                        
        compositenodemacros.hh                  
        filteredcompositenode.hh                
        filters.hh                              
        fixedcapacitystack.hh                   
        generictransformationdescriptors.hh     
        leafnode.hh                             
        nodeinterface.hh                        
        nodetags.hh                             
        pairtraversal.hh                        
        powercompositenodetransformationtemplates.hh 
        powernode.hh                            
        proxynode.hh                            
        simpletransformationdescriptors.hh      
        transformation.hh                       
        transformationutilities.hh              
        traversal.hh                            
        traversalutilities.hh                   
        treepath.hh                             
        utility.hh                              
        variadiccompositenode.hh                
        visitor.hh)

set(noinst_LTLIBRARIES  libpdelabtypetree.la)
add_library(pdelabtypetree OBJECT
        utility.cc)


# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${typetree_HEADERS} DESTINATION ${typetreedir})