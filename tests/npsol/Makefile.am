#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = tstillum tstrosen tstlinear tstnonlinear
check_PROGRAMS = $(TESTS)

tstillum_SOURCES = tstillum.C illum.C MAIN.c tstfcn.h
tstrosen_SOURCES = tstrosen.C rosen.C MAIN.c tstfcn.h
tstlinear_SOURCES = tstlinear.C rosen.C MAIN.c tstfcn.h
tstnonlinear_SOURCES = tstnonlinear.C sixtyfive.C MAIN.c tstfcn.h

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Provide libraries to be linked in.

tstillum_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstrosen_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstlinear_LDADD = $(top_builddir)/lib/libopt.la \
		  $(top_builddir)/lib/libnewmat.la \
		  $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
tstnonlinear_LDADD = $(top_builddir)/lib/libopt.la \
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
