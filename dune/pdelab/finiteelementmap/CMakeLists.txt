set(finiteelementmapdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/finiteelementmap)

set(finiteelementmap_HEADERS
        bdm1cube2dfem.hh
        bdm1simplex2dfem.hh
        brezzidouglasmarinifem.hh
        edges0.5fem.hh
        finiteelementmap.hh
        global.hh
        l2orthonormal.hh
        mimeticfem.hh
        monomfem.hh
        opbfem.hh
        p0fem.hh
        pk1dbasis.hh
        pkfem.hh
        pkqkfem.hh
        powerfem.hh
        qkdg.hh
        qkdggl.hh
        qkfem.hh
        rannacher_turek2dfem.hh
        rannacherturekfem.hh
        raviartthomasfem.hh
        rt0cube2dfem.hh
        rt0cube3dfem.hh
        rt0simplex2dfem.hh
        rt1cube2dfem.hh
        rt1cube3dfem.hh
        rt1simplex2dfem.hh
        rt2cube2dfem.hh
        variablemonomfem.hh
        variableopbfem.hh
        variableqkdgfem.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${finiteelementmap_HEADERS} DESTINATION ${finiteelementmapdir})