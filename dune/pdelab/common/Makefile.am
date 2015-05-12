commondir = $(includedir)/dune/pdelab/common
common_HEADERS =				\
	benchmarkhelper.hh                      \
	borderindexidcache.hh			\
	clock.hh				\
	crossproduct.hh				\
	dofindex.hh				\
	elementmapper.hh			\
	exceptions.hh				\
	function.hh				\
	functionutilities.hh			\
	functionwrappers.hh			\
	geometrywrapper.hh			\
	globaldofindex.hh			\
	hostname.hh				\
	jacobiantocurl.hh			\
	logtag.hh				\
	multiindex.hh				\
	partitioninfoprovider.hh		\
	polymorphicbufferwrapper.hh		\
	range.hh				\
	simpledofindex.hh			\
	topologyutility.hh			\
	typetraits.hh				\
	utility.hh				\
	vtkexport.hh			\
	instationaryfilenamehelper.hh

noinst_LTLIBRARIES = libpdelabcommon.la
libpdelabcommon_la_SOURCES =			\
	clock.cc				\
	hostname.cc				\
	logtag.cc
libpdelabcommon_la_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(POSIX_CLOCK_CPPFLAGS)
libpdelabcommon_la_LDFLAGS = $(AM_LDFLAGS)	\
	$(POSIX_CLOCK_LDFLAGS)
libpdelabcommon_la_LIBADD =			\
	$(DUNE_COMMON_LIBS)			\
	$(POSIX_CLOCK_LIBS)

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
