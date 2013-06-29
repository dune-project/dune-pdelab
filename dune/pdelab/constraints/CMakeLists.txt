set(constraintsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/constraints)
set(constraints_HEADERS
    conformingconstraints.hh
    hangingnodeconstraints.hh
    hangingnodemanager.hh
    interiornodeconstraints.hh
    mdfconstraints.hh
    p0constraints.hh
    p0ghostconstraints.hh
    rt0constraints.hh)

set(SUBDIRS
    common)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

install(FILES ${constraints_HEADERS} DESTINATION ${constraintsdir})