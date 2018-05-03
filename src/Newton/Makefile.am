#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set convenience library name and include source files in this
# directory.

noinst_LTLIBRARIES = libnewton.la
libnewton_la_SOURCES = OptBaNewton.C		OptBaQNewton.C	   \
		       OptBCEllipsoid.C		OptBCFDNewton.C	   \
		       OptBCNewton.C		OptBCNewtonLike.C  \
		       OptBCQNewton.C		OptCG.C		   \
		       OptConstrFDNewton.C	OptConstrNewton.C  \
		       OptConstrNewtonLike.C	OptConstrQNewton.C \
		       OptDHNIPS.C		OptFDNewton.C	   \
		       OptFDNIPS.C		OptLBFGS.C	   \
		       OptNewton.C		OptNewtonLike.C	   \
		       OptNIPS.C		OptNIPSLike.C	   \
		       OptQNewton.C		OptQNIPS.C
if HAVE_NPSOL
libnewton_la_SOURCES += OptNPSOL.C npsol_setup.c
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
