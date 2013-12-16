set(adaptivitydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/adaptivity)
set(adaptivity_HEADERS  adaptivity.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${adaptivity_HEADERS} DESTINATION ${adaptivitydir})