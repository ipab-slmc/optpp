#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Set convenience library name and include source files in this
# directory.

bin_PROGRAMS = optpp
optpp_SOURCES = AppLauncher.C	NewtonProblem.C	NPSOLProblem.C	\
		Problem.C	CGProblem.C	NIPSProblem.C	\
		opt++.C		VariableList.C	FDNLF1APP.C	\
		NLF0APP.C	PDSProblem.C

# Provide location of additional include files.

INCLUDES = -I$(top_srcdir)/newmat11 -I$(top_srcdir)/include


# Provide libraries to be linked in.

optpp_LDADD = $(top_srcdir)/lib/libopt.la \
	      $(top_srcdir)/lib/libnewmat.la \
	      $(NPSOL_LIB) $(BLAS_LIBS) $(FLIBS)

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL \
	     Opt.xml ReadOptAction.java SaveOptAction.java

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
