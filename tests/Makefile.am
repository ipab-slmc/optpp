#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Source subdirectories to be included during the build in order to
# perform regression tests.

SUBDIRS = constraints hock parallel uncon
if HAVE_XML
SUBDIRS += xml
endif
if HAVE_NPSOL
SUBDIRS += npsol
endif
DIST_SUBDIRS = constraints hock parallel uncon xml npsol

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in
