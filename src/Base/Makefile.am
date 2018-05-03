#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set convenience library name and include source files in this
# directory.

noinst_LTLIBRARIES = libbase.la
libbase_la_SOURCES = Appl_Data.C	Appl_Data_NPSOL.C \
		     backtrack.C	dogleg.C	  \
		     FDNLF1.C		linesearch.C	  \
		     LSQNLF.C		mcsrch.C	  \
		     NLF0.C		NLF1.C		  \
		     NLF2.C		NLP0.C		  \
		     NLP1.C		NLP2.C		  \
		     NLP.C		TOLS.C		  \
		     trustpds.C		trustregion.C

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
