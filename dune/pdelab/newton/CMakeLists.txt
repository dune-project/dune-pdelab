set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/newton)
set(common_HEADERS  newton.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})