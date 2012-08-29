set(constraintsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/constraints)
set(constraints_HEADERS                    
        constraints.hh                          
        constraintsparameters.hh                        
        constraintstransformation.hh    
        interiornodeconstraints.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${constraints_HEADERS} DESTINATION ${constraintsdir})