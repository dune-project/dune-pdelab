add_subdirectory(typetree)
set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/common)
set(common_HEADERS
        benchmarkhelper.hh
        clock.hh                                
        countingptr.hh                          
        cpstoragepolicy.hh                      
        crossproduct.hh                         
        exceptions.hh
        function.hh                             
        functionutilities.hh                    
        functionwrappers.hh                     
        geometrywrapper.hh                      
        hostname.hh                             
        jacobiantocurl.hh                       
        logtag.hh                               
        multiindex.hh                           
        range.hh
	topologyutility.hh
        typetree.hh                             
        typetraits.hh                           
        vtkexport.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules


dune_add_library(dunepdelab
        clock.cc                                
        hostname.cc                             
        logtag.cc
	_DUNE_TARGET_OBJECTS:pdelabtypetree_
	ADD_LIBS ${DUNE_LIBS}
	)

install(FILES ${common_HEADERS} DESTINATION ${commondir})
