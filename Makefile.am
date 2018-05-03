#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

ACLOCAL_AMFLAGS = -I config

# Source subdirectories to be included during the build.

SUBDIRS = newmat11 src tests
if HAVE_DOXYGEN
SUBDIRS += docs
endif
DIST_SUBDIRS = newmat11 src tests docs

# Header files to be included in the distribution.

include_HEADERS = include/abbrev_copyright.h	include/AppLauncher.h	     \
		  include/Appl_Data.h		include/Appl_Data_NPSOL.h    \
		  include/BoolVector.h		include/BoundConstraint.h    \
		  include/cblas.h		include/CGProblem.h	     \
		  include/common.h		include/CompoundConstraint.h \
		  include/ConstraintBase.h	include/Constraint.h	     \
		  include/GenSetBase.h		include/GenSetBox2d.h	     \
		  include/GenSet.h		include/GenSetMin.h	     \
		  include/GenSetStd.h		include/globals.h	     \
		  include/ioformat.h		include/LinearConstraint.h   \
		  include/LinearEquation.h	include/LinearInequality.h   \
		  include/LSQNLF.h		include/NewtonProblem.h	     \
		  include/NIPSProblem.h		include/NLFAPP.h	     \
		  include/NLF.h			include/NLP0.h		     \
		  include/NLP1.h		include/NLP2.h		     \
		  include/NLPBase.h		include/NLP.h		     \
		  include/NonLinearConstraint.h include/NonLinearEquation.h  \
		  include/NonLinearInequality.h include/NPSOLProblem.h	     \
		  include/OptBaNewton.h		include/OptBaQNewton.h	     \
		  include/OptBCEllipsoid.h	include/OptBCFDNewton.h	     \
		  include/OptBCNewton.h		include/OptBCNewtonLike.h    \
		  include/OptBCQNewton.h	include/OptCG.h		     \
		  include/OptConstrFDNewton.h	include/OptConstrNewton.h    \
		  include/OptConstrNewtonLike.h	include/OptConstrQNewton.h   \
		  include/OptDHNIPS.h		include/OptDirect.h	     \
		  include/OptFDNewton.h		include/OptFDNIPS.h	     \
		  include/OptGSS.h		include/Opt.h		     \
		  include/OptLBFGS.h		include/OptNewton.h	     \
		  include/OptNewtonLike.h	include/OptNIPS.h	     \
		  include/OptNIPSLike.h		include/OptNPSOL.h	     \
		  include/Opt_PARAMS.h		include/OptPDS.h	     \
		  include/OptppArray.h		include/OptppExceptions.h    \
		  include/OptppFatalError.h	include/OptppSmartPtr.h	     \
		  include/OptQNewton.h		include/OptQNIPS.h	     \
		  include/pds.h			include/PDSProblem.h	     \
		  include/Problem.h		include/proto.h		     \
		  include/TOLS.h		include/VariableList.h

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL INSTALL README  \
	     include/COPYRIGHT include/GNU_LGPL \
	     lib/optMaui.jar \
             docs/COPYRIGHT  docs/GNU_LGPL  \
             docs/images docs/techreport docs/guide    \
             docs/Makefile.am  docs/optpp_doc.conf    


# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *~ bin/*~ config/*~ include/*~ lib/*~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = aclocal.m4 configure Makefile.in	       \
		       config/config.guess config/config.sub   \
		       config/depcomp config/install-sh	       \
		       config/ltmain.sh config/missing	       \
		       include/OPT++_config.h.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -fr ii_files SunWS_cache so_locations
