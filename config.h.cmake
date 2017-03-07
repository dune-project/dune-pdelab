/* begin dune-pdelab
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/
/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of dune-pdelab */
#define DUNE_PDELAB_VERSION "${DUNE_PDELAB_VERSION}"

/* Define to the major version of dune-pdelab */
#define DUNE_PDELAB_VERSION_MAJOR ${DUNE_PDELAB_VERSION_MAJOR}

/* Define to the minor version of dune-pdelab */
#define DUNE_PDELAB_VERSION_MINOR ${DUNE_PDELAB_VERSION_MINOR}

/* Define to the revision of dune-pdelab */
#define DUNE_PDELAB_VERSION_REVISION ${DUNE_PDELAB_VERSION_REVISION}

/* This is only true if PETSc was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#ifndef HAVE_PETSC
#cmakedefine HAVE_PETSC ENABLE_PETSC
#endif

/* This is only true if Eigen3 was found by configure */
#ifndef HAVE_EIGEN
#cmakedefine HAVE_EIGEN ENABLE_EIGEN
#endif

/* Define to 1 if GCC's __typeof__ extension is supported. */
#ifndef HAVE_GCC___TYPEOF__
#cmakedefine HAVE_GCC___TYPEOF__ 1
#endif

/* Define to 1 if sequential UG has been found */
#ifndef PDELAB_SEQUENTIAL_UG
#cmakedefine PDELAB_SEQUENTIAL_UG 1
#endif

/* end dune-pdelab */
