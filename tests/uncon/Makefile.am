#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tstqnewton tstgnewton tstnewton tstfdnlf1 tstcg tstLBFGS
check_PROGRAMS = $(TESTS)

tstqnewton_SOURCES = tstqnewton.C rosen.C tstfcn.h
tstgnewton_SOURCES = tstgnewton.C rosen.C tstfcn.h
tstnewton_SOURCES = tstnewton.C rosen.C tstfcn.h
tstfdnlf1_SOURCES = tstfdnlf1.C rosen.C tstfcn.h
tstcg_SOURCES = tstcg.C rosen.C tstfcn.h
tstLBFGS_SOURCES = tstLBFGS.C rosen.C tstfcn.h

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Provide libraries to be linked in.

tstqnewton_LDADD = $(top_builddir)/lib/libopt.la \
		   $(top_builddir)/lib/libnewmat.la \
		   $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstgnewton_LDADD = $(top_builddir)/lib/libopt.la \
		   $(top_builddir)/lib/libnewmat.la \
		   $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstnewton_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstfdnlf1_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstcg_LDADD = $(top_builddir)/lib/libopt.la \
	      $(top_builddir)/lib/libnewmat.la \
	      $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstLBFGS_LDADD = $(top_builddir)/lib/libopt.la \
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
