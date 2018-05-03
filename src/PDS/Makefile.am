#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set convenience library name and include source files in this
# directory.

noinst_LTLIBRARIES = libpds.la
libpds_la_SOURCES = create_scheme.C	dqrdc.c	  \
		    make_search.c	pds.C	  \
		    pdschk.C		pdscld.c  \
		    pdscom.c		pdsdgn.c  \
		    pdsdone.c		pdseql.c  \
		    pdsget.c		pdshrk.c  \
		    pdsinit.C		pdslen.c  \
		    pdsopt.C		pdsrgt.c  \
		    pdsstep.C		pdsupd.c  \
		    pdswork.C		sort.c	  \
		    writes.c
if HAVE_MPI
libpds_la_SOURCES += pdsglb.c pdsgop.c pdswap.c
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
