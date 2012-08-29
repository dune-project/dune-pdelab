# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
set(NORMALTESTS  )
set(MOSTLYCLEANFILES )
set(check_SCRIPTS  )

# programs just to build when "make check" is used
set(check_PROGRAMS  ${NORMALTESTS})

# UG flags must be in front of ALUGrid flags, else when
# dune/grid/uggrid/ugincludes.hh does #include "parallel.h",
# ALUGrid-1.14-inst/include/serial/parallel.h will be included...
set(AM_CPPFLAGS  ${AM_CPPFLAGS}                             
        ${UG_CPPFLAGS}                                  
        ${ALUGRID_CPPFLAGS}                             
        ${AMIRAMESH_CPPFLAGS}                           
        ${DUNE_CPPFLAGS}                                
        ${DUNEMPICPPFLAGS}                              
        ${PARMETIS_CPPFLAGS})
set(AM_LDFLAGS  ${AM_LDFLAGS}                       
        ${UG_LDFLAGS}                           
        ${ALUGRID_LDFLAGS}                      
        ${AMIRAMESH_LDFLAGS}                    
        ${DUNE_LDFLAGS}                         
        ${DUNEMPILDFLAGS}                       
        ${PARMETIS_LDFLAGS})
set(LDADD                                          
        ${PARMETIS_LIBS}                        
        ${ALUGRID_LIBS}                         
        ${AMIRAMESH_LIBS}                       
        ${DUNE_LIBS}                            
        ${DUNEMPILIBS}                          
        ${UG_LIBS}                              
        ${LDADD})

#noinst_HEADERS =

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# TRUNK TARGETS FOLLOW

list(APPEND NORMALTESTS poisson-old)
add_executable(poisson_old poisson-old.cc)
list(APPEND MOSTLYCLEANFILES poisson-old_*.vtu poisson-old_*timings.txt)

list(APPEND NORMALTESTS laplacedirichletccfv-old)
add_executable(laplacedirichletccfv_old laplacedirichletccfv-old.cc)
list(APPEND MOSTLYCLEANFILES laplacedirichletccfv-old_*.vtu laplacedirichletccfv-old_*timings.txt)

foreach(_test ${NORMALTEST})
  add_dune_ug_flags(${_test})
  add_dune_alugrid_flags(${_test})
  add_dune_amiramesh_flags(${_test})
  add_dune_mpi_flags(${_test})
  add_dune_parmetis_flags(${_test})
endforeach(_test ${NORMALTEST})

# THIS IS A SEPARATION MARKER TO MINIMIZE SVN MERGE CONFLICTS
# BRANCH TARGETS FOLLOW

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

set(DISTCHECK_CONFIGURE_FLAGS  --with-dune=${DUNEROOT} CXX="${CXX}" CC="${CC}")

# list of tests to run
set(TESTS  ${NORMALTESTS})

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})


# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})