/* begin dune-pdelab
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/
/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME"

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

/* This is only true if PETSc was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#cmakedefine HAVE_PETSC ENABLE_PETSC

/* Define to 1 if std::initializer_list is supported. */
#cmakedefine HAVE_INITIALIZER_LIST 1

/* Define to 1 if you have the <tr1/unordered_set> header file. */
#cmakedefine HAVE_TR1_UNORDERED_SET 1

/* Define to 1 if you have the <tr1/unordered_map> header file. */
#cmakedefine HAVE_TR1_UNORDERED_MAP 1

/* Define to 1 if you have the <unordered_map> header file. */
#cmakedefine HAVE_UNORDERED_MAP 1

/* Define to 1 if you have the <unordered_set> header file. */
#cmakedefine HAVE_UNORDERED_SET 1

/* Define to 1 if decltype if supported. */
#cmakedefine HAVE_STD_DECLTYPE 1

/* Define to 1 if GCC's __typeof__ extension is supported. */
#cmakedefine HAVE_GCC___TYPEOF__ 1

/* end dune-pdelab */
