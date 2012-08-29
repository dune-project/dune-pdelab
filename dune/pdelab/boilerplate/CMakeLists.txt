set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/boilerplate)
set(common_HEADERS  pdelab.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})