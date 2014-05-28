set(M4FILES
  dune-pdelab.m4
  dune-posix-clock.m4
  eigen.m4
  petsc.m4)

set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/dune/aclocal)
set(aclocal_DATA  ${M4FILES})

set(EXTRA_DIST  ${M4FILES})

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})
