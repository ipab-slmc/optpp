dnl OPT++ Options

AC_DEFUN([OPTPP_OPTIONS],[

dnl Check for BLAS library.

   AC_F77_LIBRARY_LDFLAGS
   ACX_BLAS
   AM_CONDITIONAL([HAVE_BLAS], [test "x$acx_blas_ok" = xyes])

dnl Check for and set up MPI to build parallel OPT++.

   have_mpi=no
   AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi],
			             [build parallel version of OPT++]),
		[enable_mpi=$enableval], [enable_mpi=no])

   if test "x$enable_mpi" = xyes; then
      ACX_MPI([CXX="$MPICXX" AC_LANG_PUSH([C]) ACX_MPI([CC="$MPICC"
	       LIBS="$MPILIBS $LIBS" have_mpi=yes]) AC_LANG_POP([C])])

    dnl GM-MPI option check.
#    AC_ARG_ENABLE([gm],AS_HELP_STRING([--enable-gm],[turn GM support on]),
#		  [enable_gm=$enableval],[enable_gm=no])
#    if test "x$enable_gm" = xyes; then
#      AC_CHECK_LIB([gm],[gm_init],MPILIBS="$MPILIBS -lgm")
#    fi
   fi

   if test "x$have_mpi" = xyes; then
      AC_CHECK_LIB([mpich], [MPI_Get_version], [CXXFLAGS="$CXXFLAGS -DMPICH_IGNORE_CXX_SEEK"])
      AC_DEFINE(WITH_MPI, 1, [Define if you are building parallel OPT++.])
      AC_DEFINE(SHARED, 1, [Define if you have a shared file system.])
   fi
   AM_CONDITIONAL([HAVE_MPI], [test "x$have_mpi" = xyes])

   have_xml=no
   AM_CONDITIONAL([HAVE_XML], [test "x$have_xml" = xyes])

dnl Check for and set up to build with NPSOL.

   have_npsol=no
   AC_ARG_WITH(npsol, [AC_HELP_STRING([--with-npsol=<lib>],
				      [use NPSOL library <lib>])])

   case $with_npsol in
      yes) ;;
      no | "") have_npsol=disable ;;
      -* | */* | *.a | *.so | *.so.* | *.o) NPSOL_LIB="$with_npsol" ;;
      *) NPSOL_LIB="-l$with_npsol" ;;
   esac

   if test $have_npsol = no; then
      if test "x$NPSOL_LIB" = x; then
         NPSOL_LIB=""
      else
         AC_MSG_CHECKING([for $NPSOL_LIB])
	 if test -f $NPSOL_LIB; then
            have_npsol=yes
	 else
	    NPSOL=""
	 fi
         AC_MSG_RESULT($have_npsol)
      fi
   fi

   if test $have_npsol = no; then
      AC_MSG_CHECKING([for /usr/local/lib/libnpsol.a])
      if test -f "/usr/local/lib/libnpsol.a"; then
	 have_npsol=yes
	 NPSOL_LIB="/usr/local/lib/libnpsol.a"
      else
	 NPSOL_LIB=""
      fi
      AC_MSG_RESULT($have_npsol)
   fi

   if test $have_npsol = no; then
      AC_MSG_CHECKING([for /usr/local/lib/npsol.a])
      if test -f "/usr/local/lib/npsol.a"; then
	 have_npsol=yes
	 NPSOL_LIB="/usr/local/lib/npsol.a"
      else
	 NPSOL_LIB=""
      fi
      AC_MSG_RESULT($have_npsol)
   fi

   if test $have_npsol = no; then
      AC_MSG_CHECKING([for /usr/lib/libnpsol.a])
      if test -f "/usr/lib/libnpsol.a"; then
	 have_npsol=yes
	 NPSOL_LIB="/usr/lib/libnpsol.a"
      else
	 NPSOL_LIB=""
      fi
      AC_MSG_RESULT($have_npsol)
   fi

   if test $have_npsol = no; then
      AC_MSG_CHECKING([for /usr/lib/npsol.a])
      if test -f "/usr/lib/npsol.a"; then
	 have_npsol=yes
	 NPSOL_LIB="/usr/lib/npsol.a"
      else
	 NPSOL_LIB=""
      fi
      AC_MSG_RESULT($have_npsol)
   fi

   AC_SUBST(NPSOL_LIB)
   AM_CONDITIONAL([HAVE_NPSOL], [test "x$have_npsol" = xyes])

dnl Check for doxygen to build HTML documentation.

   AC_ARG_ENABLE(html-docs, AC_HELP_STRING([--enable-html-docs],
			    [build HTML documentation using doxygen]),
		[enable_html_docs=$enableval], [enable_html_docs=no])

   if test "x$enable_html_docs" = xyes; then
      AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)
      if test "x$DOXYGEN" = x; then
	 have_doxygen=no
      else
	 have_doxygen=yes
      fi
   else
      have_doxygen=no
   fi

   AM_CONDITIONAL([HAVE_DOXYGEN], [test "x$have_doxygen" = xyes])

])
