#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

# Header files to be included in the distribution.

include_HEADERS = controlw.h	include.h	myexcept.h \
		  newmatap.h	newmat.h	newmatio.h \
		  newmatnl.h	newmatrc.h	newmatrm.h \
		  precisio.h	solution.h	tmt.h

# Set main library name, include source files in this directory,
# and combine convenience libraries into main library.

lib_LTLIBRARIES = @top_builddir@/lib/libnewmat.la
@top_builddir@_lib_libnewmat_la_SOURCES = fft.C		evalue.C   \
					submat.C	cholesky.C \
					hholder.C	sort.C     \
					newmatrm.C	jacobi.C   \
					svd.C     	newmat8.C  \
					newmat7.C 	newmat6.C  \
					newmat5.C 	newmat3.C  \
					newmat4.C 	newmat2.C  \
					newmat1.C 	bandmat.C  \
					myexcept.C	newmatex.C \
					newfft.C

# Additional files to be included in the distribution.

EXTRA_DIST = newmat9.C newmatnl.C solution.C example.cpp garch.cpp   \
	     nl_ex.cpp nm_ex1.cpp nm_ex2.cpp sl_ex.cpp test_exc.cpp  \
	     tmt1.cpp tmt2.cpp tmt3.cpp tmt4.cpp tmt5.cpp tmt6.cpp   \
	     tmt7.cpp tmt8.cpp tmt9.cpp tmt1.cpp tmta.cpp tmtb.cpp   \
	     tmtc.cpp tmt.cpp tmtd.cpp tmte.cpp tmtf.cpp tmtg.cpp    \
	     tmth.cpp tmti.cpp tmtj.cpp tmtk.cpp tmtl.cpp tmtm.cpp   \
	     nm_b55.mak nm_b56.mak nm_cc.mak nm_gnu.mak nm_i5.mak    \
	     nm_il5.mak nm_m6.mak nm_ow.mak example.txt garch.txt    \
	     nl_ex.txt nm_ex1.txt nm_ex2.txt nm_targ.txt readme.txt  \
	     sl_ex.txt test_exc.txt tmt.txt rbd.css configure.in     \
	     nm11.htm CONDITIONS_OF_USE

# Files to remove by 'make distclean'

DISTCLEANFILES = *.log *.out *.ti *~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -rf ii_files SunWS_cache so_locations
