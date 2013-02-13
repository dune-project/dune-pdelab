#the dune-pdelab library
lib_LTLIBRARIES = libdunepdelab.la

# this is just a renaming of libpdelabcommon to prevent name clashes
libdunepdelab_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libdunepdelab_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libdunepdelab_la_LIBADD = ../dune/pdelab/common/libpdelabcommon.la

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
