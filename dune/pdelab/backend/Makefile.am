commondir = $(includedir)/dune/pdelab/backend

SUBDIRS = common eigen simple istl

common_HEADERS = backendselector.hh             \
		 eigen.hh			\
                 eigenmatrixbackend.hh          \
                 eigensolverbackend.hh          \
                 eigenvectorbackend.hh          \
                 istl.hh                        \
                 istlmatrixbackend.hh           \
                 istlsolverbackend.hh           \
                 istlvectorbackend.hh           \
                 novlpistlsolverbackend.hh      \
                 ovlpistlsolverbackend.hh       \
                 seqistlsolverbackend.hh        \
                 simple.hh                      \
                 solver.hh                      \
                 tags.hh

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
