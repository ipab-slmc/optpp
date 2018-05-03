#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tstfdnewtpds tstnewtpds tstpds tsttrpds tstGSS
check_PROGRAMS = $(TESTS)
#check_SCRIPTS  = processTestResults.sh

tstfdnewtpds_SOURCES = tstfdnewtpds.C SetupTest.C tstfcn.C tstfcn.h
tstnewtpds_SOURCES = tstnewtpds.C SetupTest.C tstfcn.C tstfcn.h
tstpds_SOURCES = tstpds.C SetupTest.C tstfcn.C tstfcn.h
tsttrpds_SOURCES = tsttrpds.C SetupTest.C tstfcn.C tstfcn.h
tstGSS_SOURCES = tstGSS.C SetupTest.C tstfcn.C tstfcn.h

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Provide libraries to be linked in.

tstfdnewtpds_LDADD = $(top_builddir)/lib/libopt.la \
		     $(top_builddir)/lib/libnewmat.la \
		     $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstnewtpds_LDADD = $(top_builddir)/lib/libopt.la \
		   $(top_builddir)/lib/libnewmat.la \
		   $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstpds_LDADD = $(top_builddir)/lib/libopt.la \
	       $(top_builddir)/lib/libnewmat.la \
	       $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tsttrpds_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstGSS_LDADD = $(top_builddir)/lib/libopt.la \
	       $(top_builddir)/lib/libnewmat.la \
	       $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = myscheme *.log *.out* *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
