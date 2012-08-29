set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/linearproblem)

set(my_HEADERS                                     
        stationarymatrix.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${my_HEADERS} DESTINATION ${mydir})