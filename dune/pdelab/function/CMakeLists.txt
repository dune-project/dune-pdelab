set(myincdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/function)

set(myinc_HEADERS
        localfunction.hh
        localfunctionhelper.hh
        oldinterfaceadapter.hh
        tags.hh
        # olf interface ...
        const.hh
        division.hh
        inverse.hh
        memberadaptor.hh
        minus.hh
        product.hh
        scalarscaled.hh
        selectcomponent.hh
        sqr.hh
        sqrt.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules


install(FILES ${myinc_HEADERS} DESTINATION ${myincdir})