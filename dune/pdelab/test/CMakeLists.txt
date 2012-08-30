# $Id: Makefile.am 586 2006-04-28 13:41:31Z mblatt $

#defined empty so we can add to it later
set(NORMALTESTS  )
set(MOSTLYCLEANFILES )
set(check_SCRIPTS  )

# programs just to build when "make check" is used
set(check_PROGRAMS  ${NORMALTESTS})

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
target_link_libraries(testclock dunepdelab ${DUNE_LIBS})

#add compile and link flags
# UG flags must be in front of ALUGrid flags, else when
# dune/grid/uggrid/ugincludes.hh does #include "parallel.h",
# ALUGrid-1.14-inst/include/serial/parallel.h will be included...
foreach(_test ${NORMALTEST})
  add_dune_ug_flags(${_test})
  add_dune_alugrid_flags(${_test})
  add_dune_amiramesh_flags(${_test})
  add_dune_mpi_flags(${_test})
  add_dune_parmetis_flags(${_test})
endforeach(_test ${NORMALTEST})
# don't include all the grid stuff

# The following lines are comment out as the translation is not clear. TODO review
#
# set(testclock_ LDFLAGS ${AM_LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(testclock_ LDADD ${LDADD})


# list of tests to run
set(TESTS  ${NORMALTESTS})


foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})


# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})