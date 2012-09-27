set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/adaptivity)
set(common_HEADERS  adapt.hh adaptivity.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})