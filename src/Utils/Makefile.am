#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set convenience library name and include source files in this
# directory.

noinst_LTLIBRARIES = libutils.la
libutils_la_SOURCES = BoolVector.C		file_cutils.c	  \
		      ioformat.C		mcholesky.C	  \
		      OptppExceptions.C		OptppFatalError.C \
		      print.C			timers.c
if !HAVE_BLAS
libutils_la_SOURCES += linalg.c
endif

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
