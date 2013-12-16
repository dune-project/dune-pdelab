set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/finiteelement)
set(my_HEADERS
        localbasiscache.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${my_HEADERS} DESTINATION ${mydir})