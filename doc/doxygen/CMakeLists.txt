# $Id: Makefile.am 4161 2008-05-20 21:47:51Z christi $

# setting like in dune-web
set(BASEDIR ../..)
set(CURDIR doc/doxygen)
set(WHAT )

set(doxygen_DATA  ${WHAT})

#we need this in the distribution
set(EXTRA_DIST  ${WHAT} ${DOXYGENDISTFILES})

# include not needed for CMake
# include $(top_srcdir)/am/doxygen
# include not needed for CMake
# include $(top_srcdir)/am/global-rules

# file to clean only in svn and not in tarball tree
set(SVNCLEANFILES  html ${WHAT} ${DOXYGENHEADER})


# Create and install doxygen documentation
add_doxygen_target()