dnl -*- mode: autoconf; tab-width: 8; indent-tabs-mode: nil; -*-
dnl vi: set et ts=8 sw=2 sts=2:
# $Id$
# searches for PETSc headers and libs

AC_DEFUN([DUNE_PATH_PETSC],[
  AC_MSG_CHECKING(for PETSc library)
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])
  AC_REQUIRE([DUNE_MPI])

  #
  # USer hints ...
  #
  AC_LANG_PUSH([C++])
  AC_ARG_VAR([PETSC], [PETSc library location])
  AC_ARG_WITH([petsc],
    [AC_HELP_STRING([--with-petsc],[user defined path to PETSc library])],
    [
      # --with-petsc supersedes $PETSC
      PETSC=""
      if test "$withval" != no ; then
          if test -d "$withval" ; then
	    # get absolute path
	    with_petsc=`eval cd $withval 2>&1 && pwd`
            include_path=include
            lib_path=lib
            echo "$with_petsc/$include_path/petsc.h"
            if test -f "$with_petsc/$include_path/petsc.h" ; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
            fi
          else
            with_petsc="no"
	    AC_MSG_RESULT(no)
      fi
	else
	    AC_MSG_RESULT(no)
	fi
	],
    [
	if test -n "$PETSC" ; then
          if test -d "$PETSC" ; then
	    # get absolute path
	    with_petsc=`eval cd $PETSC 2>&1 && pwd`
            PETSC=""
            include_path=include
            lib_path=lib
            if test ! -f "$with_petsc/$include_path/petsc.h" ; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
            fi
          else
            PETSC=""
            with_petsc=no
            AC_MSG_RESULT(no)
          fi
        else
            with_petsc=/usr/
            include_path=include
            lib_path=lib
            if test ! -f "$with_petsc/$include_path/petsc.h" ; then
                with_petsc=/usr/local/
                if test ! -f "$with_petsc/$include_path/petsc.h" ; then
                    with_petsc="no"
                    AC_MSG_RESULT(failed)
                else
                    AC_MSG_RESULT(yes)
                fi
            else
                AC_MSG_RESULT(yes)
            fi
        fi
        ])
  

  # store old values
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  
  ## do nothing if --without-petsc is used
  if test x"$with_petsc" != x"no" ; then
          
      # defaultpath
      PETSC_LIB_PATH="$with_petsc$lib_path"
      PETSC_INCLUDE_PATH="$with_petsc/$include_path"
                  
      PETSC_LIBS="-lpetsc $LAPACK_LIBS $BLAS_LIBS -lX11 $DUNEMPILIBS"
      PETSC_LDFLAGS="-L$PETSC_LIB_PATH $DUNEMPILDFLAGS"

      # set variables so that tests can use them
      CPPFLAGS="$CPPFLAGS -I$PETSC_INCLUDE_PATH $DUNEMPICPPFLAGS"

      # check for central header
      AC_CHECK_HEADER([petsc.h],[
              PETSC_CPPFLAGS="-I$PETSC_INCLUDE_PATH $DUNEMPICPPFLAGS"
	      HAVE_PETSC="1"],[
	      HAVE_PETSC="0"
	      AC_MSG_WARN([petsc.h not found in $PETSC_INCLUDE_PATH with $CPPFLAGS])]
      )

      PETSC_CPPFLAGS="${PETSC_CPPFLAGS} -DENABLE_PETSC=1"

#      AC_LANG_PUSH([C++])
      
      # if header is found check for the libs

      LIBS="-lm $LIBS -lX11 $LAPACK_LIBS $BLAS_LIBS $X_LIBS $DUNEMPILIBS $DUNEMPILDFLAGS"
      
      if test x$HAVE_PETSC = x1 ; then
	  DUNE_CHECK_LIB_EXT([$PETSC_LIB_PATH], [petsc], [PetscTrMalloc],
              [
		  PETSC_LIBS="-L$PETSC_LIB_PATH -lpetsc $X_LIBS -lX11 $LAPACK_LIBS $BLAS_LIBS -lm $DUNEMPILIBS $DUNEMPILDFLAGS"
		  LIBS="$PETSC_LIBS $ac_save_LIBS"
              ],[
		  HAVE_PETSC="0"
		  AC_MSG_WARN(libpetsc not found!)
              ])
      fi

#      AC_LANG_POP([C++])
      
      # pre-set variable for summary
      #with_petsc="no"
      
      # did it work?
      AC_MSG_CHECKING(PETSc in $with_petsc)
      if test x$HAVE_PETSC = x1 ; then
	  AC_SUBST(PETSC_LIBS, $PETSC_LIBS)
	  AC_SUBST(PETSC_LDFLAGS, $PETSC_LDFLAGS)
	  AC_SUBST(PETSC_CPPFLAGS, $PETSC_CPPFLAGS)
	  AC_DEFINE(HAVE_PETSC,ENABLE_PETSC,[Define if you have the PETSc library.
		  This is only true if the application uses the PETSC_CPPFLAGS])
	  AC_MSG_RESULT(ok)
	  
    # add to global list
	  DUNE_PKG_LIBS="$PETSC_LIBS $DUNE_PKG_LIBS"
	  DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $PETSC_LDFLAGS"
	  DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $PETSC_CPPFLAGS"
	  
    # re-set variable correctly
	  with_petsc="yes"
      else
	  with_petsc="no"
	  AC_MSG_RESULT(failed)
      fi 
      
  # end of "no --without-petsc"
  else
  	with_petsc="no"
  fi

  AC_LANG_POP([C++]) 

  # tell automake	
  AM_CONDITIONAL(PETSC, test x$HAVE_PETSC = x1)
  
  # restore variables
  LDFLAGS="$ac_save_LDFLAGS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LIBS="$ac_save_LIBS"
  
  DUNE_ADD_SUMMARY_ENTRY([PETSc],[$with_petsc])

])
