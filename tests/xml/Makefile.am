#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set list of of tests to be built run by 'make check' and provide the
# relevant source files.

TESTS = rosen hockfcns
check_PROGRAMS = $(TESTS)

rosen_SOURCES = rosen.C
hockfcns_SOURCES = hockfcns.C

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Provide libraries to be linked in.

rosen_LDADD = $(top_builddir)/lib/libopt.la \
	      $(top_builddir)/lib/libnewmat.la \
	      $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)
hockfcns_LDADD = $(top_builddir)/lib/libopt.la \
		 $(top_builddir)/lib/libnewmat.la \
		 $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL \
	     tstCG.xml tstFD.xml tstHock14.xml tstHock65.xml \
	     tstNewton.xml tstNIPS.xml tstPDS.xml tstQNewton.xml \
	     tstTRPDS.xml

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
