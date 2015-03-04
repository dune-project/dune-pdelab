set(gridoperatortbbtdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/tbb)

set(gridoperatortbb_HEADERS
        assembler.hh
        jacobianapplyengine.hh
        jacobianengine.hh
        localassembler.hh
        residualengine.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperatortbb_HEADERS} DESTINATION ${gridoperatortbbtdir})
