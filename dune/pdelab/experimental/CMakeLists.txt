set(noinst_DATA  README.txt)
set(EXTRA_DIST  README.txt)

set(SUBDIRS                                        
        common)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
