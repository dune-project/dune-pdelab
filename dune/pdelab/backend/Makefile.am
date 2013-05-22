commondir = $(includedir)/dune/pdelab/backend

SUBDIRS = common dense istl

common_HEADERS = backendselector.hh             \
		 dense.hh			\
                 eigenmatrixbackend.hh          \
                 eigensolverbackend.hh          \
                 eigenvectorbackend.hh          \
                 istlmatrixbackend.hh           \
                 istlsolverbackend.hh           \
                 istlvectorbackend.hh           \
                 novlpistlsolverbackend.hh      \
                 ovlpistlsolverbackend.hh       \
                 petscmatrixbackend.hh          \
                 petscnestedvectorbackend.hh    \
                 petscutility.hh                \
                 petscvectorbackend.hh          \
                 seqistlsolverbackend.hh        \
                 solver.hh	                \
                 tags.hh			\
                 vectorutilities.hh

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
