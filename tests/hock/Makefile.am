#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tsthock1 tsthock2 tsthock5 tsthock6 tsthock7 tsthock10	    \
	tsthock13 tsthock14 tsthock26 tsthock28 tsthock35 tsthock65 \
	tsthock77 tsthock78
check_PROGRAMS = $(TESTS)

tsthock1_SOURCES = tsthock1.C hockfcns.C hockfcns.h
tsthock2_SOURCES = tsthock2.C hockfcns.C hockfcns.h
tsthock5_SOURCES = tsthock5.C hockfcns.C hockfcns.h
tsthock6_SOURCES = tsthock6.C hockfcns.C hockfcns.h
tsthock7_SOURCES = tsthock7.C hockfcns.C hockfcns.h
tsthock10_SOURCES = tsthock10.C hockfcns.C hockfcns.h
tsthock13_SOURCES = tsthock13.C hockfcns.C hockfcns.h
tsthock14_SOURCES = tsthock14.C hockfcns.C hockfcns.h
tsthock26_SOURCES = tsthock26.C hockfcns.C hockfcns.h
tsthock28_SOURCES = tsthock28.C hockfcns.C hockfcns.h
tsthock35_SOURCES = tsthock35.C hockfcns.C hockfcns.h
tsthock65_SOURCES = tsthock65.C hockfcns.C hockfcns.h
tsthock77_SOURCES = tsthock77.C hockfcns.C hockfcns.h
tsthock78_SOURCES = tsthock78.C hockfcns.C hockfcns.h

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Provide libraries to be linked in.

tsthock1_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock2_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock5_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock6_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock7_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock10_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock13_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock14_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock26_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock28_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock35_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock65_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock77_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsthock78_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
