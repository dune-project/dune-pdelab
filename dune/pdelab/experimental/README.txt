This "experimental" directory contains code which for various reasons should
not be used by unsuspecting users but nevertheless should be kept around.  The
following guidelines apply:

 * This directory is included in the build system so people can run "make
   headercheck" to see if a header still compiles.

 * Headers in this directory are not installed, they are mentioned in
   noinst_HEADERS in the Makefile.am.

 * Headers in this directory _are_ distributed, i.e. do not mention them in
   nodist_noinst_HEADERS.  The reason is that not distributing them would make
   the Makefile.am much more complicated.

 * The directory structure below dune/pdelab/experimental should mirror the
   one below dune/pdelab.

 * When branching for a release, the whole "experimental" directory is deleted
   from the release branch.

List of reasons why a file is in experimental:

 * common/vtkexport.hh:
   untested, except for "make headercheck".
