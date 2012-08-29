# $Id: Makefile.am 3930 2007-12-05 17:17:07Z robertk $

set(SUBDIRS  doxygen)

# setting like in dune-web
set(BASEDIR ..)
set(CURDIR doc)


foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
