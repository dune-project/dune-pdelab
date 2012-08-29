set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/stationary)
set(common_HEADERS  linearproblem.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})