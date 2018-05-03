#                      -*- Automake -*-
# Process this file with automake to produce a Makefile.in.

HTML_FILES = guide/AboutDoc.html		 \
	     guide/BoundConstrainedProblems.html \
	     guide/BoundConstraints.html	 \
	     guide/Classification.html		 \
	     guide/CompoundConstraintsDoc.html	 \
	     guide/ConstrainedProblems.html	 \
	     guide/ConstraintsDoc.html		 \
	     guide/ControlParameters.html	 \
	     guide/COPYRIGHT.html		 \
	     guide/example1.html		 \
	     guide/example1_out.html		 \
	     guide/example2.html		 \
	     guide/example2_out.html		 \
	     guide/FAQS.html			 \
	     guide/footer.html			 \
	     guide/gensetGuide-format.html	 \
	     guide/GUI_XMLDoc.html		 \
	     guide/header.html			 \
	     guide/InstallDoc.html		 \
	     guide/License.html			 \
	     guide/LinearConstraints.html	 \
	     guide/MainPage.html		 \
	     guide/MethodsDoc.html		 \
	     guide/NonLinearConstraints.html	 \
	     guide/OptionsDoc.html		 \
	     guide/ParallelOptimization.html	 \
	     guide/ReferencesDoc.html		 \
	     guide/SetUp.html			 \
	     guide/subroutine1.html		 \
	     guide/subroutine1_out.html		 \
	     guide/subroutine2.html		 \
	     guide/subroutine2_out.html		 \
	     guide/tstbcnips.html		 \
	     guide/tstbcqnewton.html		 \
	     guide/tstcg.html			 \
	     guide/tstcg_out.html		 \
	     guide/tstgss.html			 \
	     guide/tstnips.html			 \
	     guide/tstnpsol.html		 \
	     guide/tstpds.html			 \
	     guide/tstqnewton.html		 \
	     guide/tsttrpds.html		 \
	     guide/UnconstrainedObject.html	 \
	     guide/UnconstrainedProblems.html	 \
	     guide/xml1.html			 \
	     guide/xml2.html			 \
	     guide/xml3.html

IMAGE_FILES = images/bckgnd.jpg		images/csmrlogo.gif	    \
	      images/doxygen.gif	images/gensetBox.gif	    \
	      images/gensetMin.gif	images/gensetStd.gif	    \
	      images/maui_advanced.png	images/maui_application.png \
	      images/maui_arrays.png	images/maui_classpath.png   \
	      images/maui_config.png	images/maui_nips.png	    \
	      images/maui_noaction.png	images/maui_opt++.png	    \
	      images/maui_startup.png	images/maui_subroutine.png  \
	      images/sandialogo.gif	images/tBird2D.gif

REPORTS = techreport/opt++.ps

# Additional files to be included in the distribution.

EXTRA_DIST = COPYRIGHT GNU_LGPL optpp_doc.conf \
	     $(HTML_FILES) $(IMAGE_FILES) $(REPORTS)

# Files to remove by 'make distclean'.

DISTCLEANFILES = *.log *.out *~ guide/*~ images/*~ techreport/*~

# Autotools-generated files to remove by 'make maintainer-clean'.

MAINTAINERCLEANFILES = Makefile.in *~

# Make the HTML documentation using doxygen.

all-local: optpp_doc.conf $(HTML_FILES) $(IMAGE_FILES)
	-$(DOXYGEN) optpp_doc.conf

# Compiler cache directories (created on some platforms) to be removed.

clean-local:
	-rm -fr html
