SUBDIRS =                                       \
	typetree

commondir = $(includedir)/dune/pdelab/common
common_HEADERS =				\
	clock.hh				\
	countingptr.hh				\
	cpstoragepolicy.hh			\
	crossproduct.hh				\
	function.hh				\
	functionutilities.hh			\
	functionwrappers.hh			\
	geometrywrapper.hh			\
	hostname.hh				\
	jacobiantocurl.hh			\
	logtag.hh				\
	multiindex.hh				\
	range.hh				\
	typetree.hh                             \
        typetraits.hh                           \
	vtkexport.hh

noinst_LTLIBRARIES = libpdelabcommon.la
libpdelabcommon_la_SOURCES =			\
	clock.cc				\
	hostname.cc				\
	logtag.cc
libpdelabcommon_la_LIBADD =			\
	-lrt

include $(top_srcdir)/am/global-rules
