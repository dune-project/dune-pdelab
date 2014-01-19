set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/multistep)
set(my_HEADERS
        cache.hh
        method.hh
        parameter.hh)

install(FILES ${my_HEADERS} DESTINATION ${mydir})
