# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
set(NORMALTESTS  )
set(MOSTLYCLEANFILES )
set(check_SCRIPTS  )

set(noinst_HEADERS                                 
        fmt.hh                                  
        gnuplotgraph.hh                         
        gridexamples.hh                         
        l2difference.hh                         
        l2norm.hh                               
        typetreetargetnodes.hh                  
        typetreetestswitch.hh                   
        typetreetestutility.hh)


set(noinst_SCRIPTS                                 
        make_pvd.sh)

set(headercheck_IGNORE                             
        typetreetargetnodes.hh                  
        typetreetestswitch.hh                   
        typetreetestutility.hh)

set(GRIDDIM  2)

# scripts are not distributed automatically
# but since this is need for make all:
set(EXTRA_DIST  make_pvd.sh)

list(APPEND NORMALTESTS testanalytic)
add_executable(testanalytic testanalytic.cc)
target_link_libraries(testanalytic dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testanalytic )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testanalytic_ LDFLAGS ${AM_LDFLAGS}            
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testanalytic_ LDADD                            
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES channel.vtu)

list(APPEND NORMALTESTS testclock)
add_executable(testclock testclock.cc)
# don't include all the grid stuff

target_link_libraries(testclock dunepdelab ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testclock_ LDFLAGS ${AM_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testclock_ LDADD ${LDADD})

if(ALUGRID_FOUND)
list(APPEND NORMALTESTS test-composed-iis-gfs)
add_executable(test_composed_iis_gfs test-composed-iis-gfs.cc)
endif(ALUGRID_FOUND )

list(APPEND NORMALTESTS testconstraints)
add_executable(testconstraints testconstraints.cc)

target_link_libraries(testconstraints dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testconstraints )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testconstraints_ LDFLAGS ${AM_LDFLAGS}         
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testconstraints_ LDADD                         
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testconstraintsp1.vtu testconstraintspowerp1.vtu)

list(APPEND NORMALTESTS testcountingptr)
add_executable(testcountingptr testcountingptr.cc)

target_link_libraries(testcountingptr dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testcountingptr )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testcountingptr_ LDFLAGS ${AM_LDFLAGS}         
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testcountingptr_ LDADD                         
#        ${ALBERTA_LIBS}                         
#        ${LDADD})

list(APPEND NORMALTESTS testfiniteelementmap)
add_executable(testfiniteelementmap testfiniteelementmap.cc)

target_link_libraries(testfiniteelementmap dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testfiniteelementmap )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testfiniteelementmap_ LDFLAGS ${AM_LDFLAGS}    
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testfiniteelementmap_ LDADD                    
#        ${ALBERTA_LIBS}                         
#        ${LDADD})

list(APPEND NORMALTESTS testfunction)
add_executable(testfunction testfunction.cc)

target_link_libraries(testfunction dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testfunction )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testfunction_ LDFLAGS ${AM_LDFLAGS}            
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testfunction_ LDADD                            
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES multi.vtu single.vtu)

list(APPEND NORMALTESTS testgridfunctionspace)
add_executable(testgridfunctionspace testgridfunctionspace.cc)
target_link_libraries(testgridfunctionspace  dunepdelab ${DUNE_LIBS})

list(APPEND NORMALTESTS testlaplacedirichletccfv)
add_executable(testlaplacedirichletccfv testlaplacedirichletccfv.cc)

target_link_libraries(testlaplacedirichletccfv dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testlaplacedirichletccfv )
add_dune_alberta_flags(testlaplacedirichletccfv )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlaplacedirichletccfv_ LDFLAGS ${AM_LDFLAGS}        
#        ${ALBERTA_LDFLAGS}                              
#        ${SUPERLU_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlaplacedirichletccfv_ LDADD                
#        ${SUPERLU_LDFLAGS} ${SUPERLU_LIBS}      
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testlaplacedirichletccfv.vtu)

list(APPEND NORMALTESTS testlaplacedirichletp12d)
add_executable(testlaplacedirichletp12d testlaplacedirichletp12d.cc)

target_link_libraries(testlaplacedirichletp12d  dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testlaplacedirichletp12d )
add_dune_alberta_flags(testlaplacedirichletp12d )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlaplacedirichletp12d_ LDFLAGS ${AM_LDFLAGS}        
#        ${ALBERTA_LDFLAGS}                              
#        ${SUPERLU_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlaplacedirichletp12d_ LDADD                
#        ${SUPERLU_LDFLAGS} ${SUPERLU_LIBS}      
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testlaplacedirichletp12d.vtu)

list(APPEND NORMALTESTS testelasticity)
add_executable(testelasticity testelasticity.cc)
target_link_libraries(testelasticity dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testelasticity )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testelasticity_ LDFLAGS ${AM_LDFLAGS}  
#        ${SUPERLU_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testelasticity_ LDADD          
#        ${SUPERLU_LDFLAGS} ${SUPERLU_LIBS}      
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testelasticity.vtu)

list(APPEND NORMALTESTS testlocalfunctionspace)
add_executable(testlocalfunctionspace testlocalfunctionspace.cc)

target_link_libraries(testlocalfunctionspace dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testlocalfunctionspace )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlocalfunctionspace_ LDFLAGS ${AM_LDFLAGS}  
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testlocalfunctionspace_ LDADD                  
#        ${ALBERTA_LIBS}                         
#        ${LDADD})

list(APPEND NORMALTESTS testmultistep)
add_executable(testmultistep testmultistep.cc)

add_dune_mpi_flags(testmultistep )
target_link_libraries(testmultistep dunepdelab ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testmultistep_ LDFLAGS ${AM_LDFLAGS}           
#        ${DUNEMPILDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testmultistep_ LDADD                           
#        ${DUNEMPILIBS}                          
#        ${LDADD})

list(APPEND MOSTLYCLEANFILES 
        s*:p*:testmultistep_yasp_P1_1d-*.vtp            
        s*:testmultistep_yasp_P1_1d-*.pvtp              
        testmultistep_yasp_P1_1d.pvd)

list(APPEND NORMALTESTS testmultistepcached)
add_executable(testmultistepcached testmultistepcached.cc)

add_dune_mpi_flags(testmultistepcached )
target_link_libraries(testmultistepcached dunepdelab ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testmultistepcached_ LDFLAGS ${AM_LDFLAGS}     
#        ${DUNEMPILDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testmultistepcached_ LDADD                     
#        ${DUNEMPILIBS}                          
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES 
        s*:p*:testmultistepcached_yasp_P1_1d-*.vtp      
        s*:testmultistepcached_yasp_P1_1d-*.pvtp        
        testmultistepcached_yasp_P1_1d.pvd)

list(APPEND NORMALTESTS testp12dinterpolation)
add_executable(testp12dinterpolation testp12dinterpolation.cc)

target_link_libraries(testp12dinterpolation dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testp12dinterpolation GRIDDIM 2)
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testp12dinterpolation_ LDFLAGS ${AM_LDFLAGS}   
#        ${ALBERTA2D_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testp12dinterpolation_ LDADD                                   
#        ${ALBERTA2D_LIBS}                       
#        ${LDADD})

list(APPEND NORMALTESTS testpk)
add_executable(testpk testpk.cc)

target_link_libraries(testpk dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testpk )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpk_ LDFLAGS ${AM_LDFLAGS}                  
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpk_ LDADD                                  
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testpk.vtu)

list(APPEND NORMALTESTS testpoisson)
add_executable(testpoisson testpoisson.cc)

target_link_libraries(testpoisson dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testpoisson )
add_dune_alberta_flags(testpoisson )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_ LDFLAGS ${AM_LDFLAGS}             
#        ${ALBERTA_LDFLAGS}                      
#        ${SUPERLU_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_ LDADD                             
#        ${SUPERLU_LDFLAGS} ${SUPERLU_LIBS}      
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES poisson_*.vtu)

list(APPEND NORMALTESTS testpoisson-globalfe)
add_executable(testpoisson_globalfe testpoisson-globalfe.cc)

target_link_libraries(testpoisson_globalfe dunepdelab ${DUNE_LIBS})
add_dune_superlu_flags(testpoisson_globalfe )
add_dune_alberta_flags(testpoisson_globalfe )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_globalfe_ LDFLAGS ${AM_LDFLAGS}    
#        ${ALBERTA_LDFLAGS}                      
#        ${SUPERLU_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_globalfe_ LDADD                    
#        ${SUPERLU_LDFLAGS} ${SUPERLU_LIBS}      
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES poisson_globalfe_*.vtu)

list(APPEND NORMALTESTS testrt0)
add_executable(testrt0 testrt0.cc)

target_link_libraries(testrt0 dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testrt0 )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testrt0_ LDFLAGS ${AM_LDFLAGS}                 
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testrt0_ LDADD                                 
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES testrt0.vtu)

list(APPEND NORMALTESTS testrt02dgridfunctionspace)
add_executable(testrt02dgridfunctionspace testrt02dgridfunctionspace.cc)

target_link_libraries(testrt02dgridfunctionspace dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testrt02dgridfunctionspace )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testrt02dgridfunctionspace_ LDFLAGS ${AM_LDFLAGS}      
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testrt02dgridfunctionspace_ LDADD              
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES rt02dgridfunctionspace-*.vtu)

list(APPEND NORMALTESTS testutilities)
add_executable(testutilities testutilities.cc)

target_link_libraries(testutilities dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testutilities )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testutilities_ LDFLAGS ${AM_LDFLAGS}           
#        ${ALBERTA_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testutilities_ LDADD                           
#        ${ALBERTA_LIBS}                         
#        ${LDADD})
list(APPEND MOSTLYCLEANFILES edger.vtu interpolated.vtu q1.vtu taylorhood.vtu)

list(APPEND TESTS ${CMAKE_CURRENT_BINARY_DIR}/testvectorwave.conf)
list(APPEND check_PROGRAMS testvectorwave)
add_executable(testvectorwave testvectorwave.cc)
target_link_libraries(testvectorwave dunepdelab ${DUNE_LIBS})

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/testvectorwave.conf
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/testvectorwave.conf ${CMAKE_CURRENT_BINARY_DIR}/testvectorwave.conf
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/testvectorwave.conf testvectorwave
  )
list(APPEND check_SCRIPTS testvectorwave.conf)

list(APPEND MOSTLYCLEANFILES 
        s*:p*:vectorwave_*.vtu                  
        s*:vectorwave_*.pvtu                    
        vectorwave_*.pvd)

list(APPEND TESTS ${CMAKE_CURRENT_BINARY_DIR}/testvectorwavecached.conf)
list(APPEND check_SCRIPTS ${CMAKE_CURRENT_BINARY_DIR}/testvectorwavecached.conf)
list(APPEND check_PROGRAMS testvectorwavecached)
add_executable(testvectorwavecached testvectorwavecached.cc)
target_link_libraries(testvectorwavecached dunepdelab ${DUNE_LIBS})
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/testvectorwavecached.conf
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/testvectorwavecached.conf ${CMAKE_CURRENT_BINARY_DIR}/testvectorwavecached.conf
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/testvectorwavecached.conf testvectorwavecached
  )
list(APPEND check_SCRIPTS testvectorwavecached.conf)

list(APPEND MOSTLYCLEANFILES 
        s*:p*:vectorwavecached_*.vtu            
        s*:vectorwavecached_*.pvtu              
        vectorwavecached_*.pvd)

list(APPEND check_PROGRAMS testtypetree testtypetree-no-sfinae testtypetree-no-variadic testtypetree-no-rvalue-refs testtypetree-legacy)
add_executable(testtypetree testtypetree.cc)
set_property(TARGET testtypetree APPEND PROPERTY COMPILE_DEFINITIONS TEST_TYPETREE)
target_link_libraries(testtypetree dunepdelab ${DUNE_LIBS})

add_executable(testtypetree_no_sfinae testtypetree.cc)
set_property(TARGET testtypetree APPEND PROPERTY COMPILE_DEFINITIONS TEST_TYPETREE_NO_SFINAE)
target_link_libraries(testtypetree_no_sfinae dunepdelab ${DUNE_LIBS})

add_executable(testtypetree_no_variadic testtypetree.cc)
set_property(TARGET testtypetree_no_variadic APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE_NO_VARIADIC)
target_link_libraries(testtypetree_no_variadic dunepdelab ${DUNE_LIBS})

add_executable(testtypetree_no_rvalue_refs testtypetree.cc)
set_property(TARGET testtypetree_no_rvalue_refs APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE_NO_RVALUE_REFERENCES)
target_link_libraries(testtypetree_no_rvalue_refs dunepdelab ${DUNE_LIBS})

add_executable(testtypetree_legacy testtypetree.cc)
set_property(TARGET testtypetree_legacy APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE_LEGACY)
target_link_libraries(testtypetree_legacy dunepdelab ${DUNE_LIBS})

list(APPEND check_PROGRAMS testtypetreetransformation)
add_executable(testtypetreetransformation testtypetreetransformation.cc)
set_property(TARGET testtypetreetransformation APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE)
target_link_libraries(testtypetreetransformation dunepdelab ${DUNE_LIBS})

list(APPEND check_PROGRAMS testinterpolate)
add_executable(testinterpolate testinterpolate.cc)
target_link_libraries(testinterpolate dunepdelab ${DUNE_LIBS})

list(APPEND check_PROGRAMS testfilteredcompositenode)
add_executable(testfilteredcompositenode testfilteredcompositenode.cc)
set_property(TARGET testfilteredcompositenode APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE)
target_link_libraries(testfilteredcompositenode dunepdelab ${DUNE_LIBS})

list(APPEND check_PROGRAMS testproxynode)
add_executable(testproxynode testproxynode.cc)
target_link_libraries(testproxynode dunepdelab ${DUNE_LIBS})
set_property(TARGET testproxynode APPEND PROPERTY
  COMPILE_DEFINITIONS TEST_TYPETREE)
target_link_libraries(testproxynode dunepdelab ${DUNE_LIBS})

if(PETSC_FOUND)
list(APPEND check_PROGRAMS testpoisson-petsc)
add_executable(testpoisson_petsc testpoisson-petsc.cc)

target_link_libraries(testpoisson_petsc dunepdelab ${DUNE_LIBS})
add_dune_petsc_flags(testpoisson_petsc )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_petsc_ LDFLAGS ${AM_LDFLAGS} ${PETSC_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testpoisson_petsc_ LDADD ${LDADD} ${PETSC_LIBS})
endif(PETSC_FOUND )

if(GMP_FOUND)
list(APPEND check_PROGRAMS testopbfem)
add_executable(testopbfem testopbfem.cc)

target_link_libraries(testopbfem  dunepdelab ${DUNE_LIBS})
add_dune_gmp_flags(testopbfem )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testopbfem_ LDFLAGS ${AM_LDFLAGS} ${GMP_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testopbfem_ LDADD ${LDADD} ${GMP_LIBS})
endif(GMP_FOUND )

list(APPEND check_PROGRAMS testrtfem)
add_executable(testrtfem testrtfem.cc)
target_link_libraries(testrtfem dunepdelab ${DUNE_LIBS})
#add_dune_alberta_flags(testrtfem)
#add_dune_ug_flags(testrtfem)
#add_dune_alugrid_flags(testrtfem)

list(APPEND check_PROGRAMS testbdmfem)
add_executable(testbdmfem testbdmfem.cc)
target_link_libraries(testbdmfem dunepdelab ${DUNE_LIBS})
add_dune_alberta_flags(testbdmfem)
add_dune_ug_flags(testbdmfem)
add_dune_alugrid_flags(testbdmfem)

list(APPEND MOSTLYCLEANFILES poisson_petsc_*.vtu)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

set(DISTCHECK_CONFIGURE_FLAGS  --with-dune=${DUNEROOT} CXX="${CXX}" CC="${CC}")


foreach(i ${TESTS} ${NORMALTESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})

# programs just to build when "make check" is used
set(check_PROGRAMS  ${NORMALTESTS})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS} ${check_PROGRAMS})