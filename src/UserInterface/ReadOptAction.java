
import Maui.Database.*;
import Maui.Interface.*;
import XML.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import javax.swing.*;


/**
 * Reads an object from storage into a maui editor.
 *
 * @version 1.0
 */

public class ReadOptAction extends MauiAction
{

 /**
	* Constructor onstructs with a label for the button, a hook to the editor, and
	* the storage system from which the object will be read.
	*
	* @see MauiAction
	*/
 public ReadOptAction()
 {
	super();
 }
	
 /**
	* Callback method, triggered when the action is invoked. 
	*/
 public void doAction(ActionEvent e, XMLObject xml)
 {

	try
	 {
		MauiFileSystem fs = null;
		
		//If we've not accessed a file through Save or Read before, then
		// use the XML to determine the directory to display.
		if(MauiFileSystem.getLatestFile().equals(""))
		 {
			String dir = config_.getAttribute("dir");
			if (dir==null)
			 {
				fs = new MauiFileSystem();
			 }
			else if (dir.equals("$HOME") || dir.equals("$(HOME)"))
			 {
				fs = new MauiFileSystem(System.getProperty("user.home"));
			 }
			else
			 {
				fs = new MauiFileSystem(dir);
			 }
		 }
		// we've got a record of the last directory and file we visited
		else 
		 {
			String latest = MauiFileSystem.getLatestFile();
			fs = new MauiFileSystem(latest);
			fs.setSelectedFile(latest);
		 }
		fs.setFilter(new ExtFileFilter("xml"));
		/* read the data into an XMLObject */
		XMLObject body = fs.retrieveObject(mainPane_);

		if (body==null)
		 return;

		XMLObject welcome = body.getChild("Welcome");
		XMLObject problem = body.getChild("ProblemSetup");
		XMLObject solver = body.getChild(2);

		String functionType = problem.getChild(1).getTag();
		String solverType = solver.getTag();

		body.addAttribute("name", "root");
		body.addAttribute("altName", "OPT");
		body.addAttribute("base", "OPT");

		welcome.addAttribute("name", "Welcome");
		welcome.addAttribute("altName", "Welcome");
		welcome.addAttribute("base", "Welcome");
		welcome.addAttribute("ownerType", "OPT");

		problem.addAttribute("name", "ProblemSetup");
		problem.addAttribute("altName", "ProblemSetup");
		problem.addAttribute("base", "ProblemSetup");
		problem.addAttribute("ownerType", "OPT");

		solver.addAttribute("name", "AlgorithmDefinition");
		solver.addAttribute("altName", solverType);
		solver.addAttribute("base", "AlgorithmDef");
		solver.addAttribute("ownerType", "OPT");

		if (functionType.equals("Application"))
		{
		    XMLObject application = problem.getChild("Application");

		    XMLObject scriptName = new XMLObject("String");
		    scriptName.addAttribute("default", application.getAttribute("scriptName"));
		    scriptName.addAttribute("name", "scriptName");
		    scriptName.addAttribute("ownerType", "Application");
		    application.removeAttribute("scriptName");
		    application.addChild(scriptName);

		    XMLObject modelDir = new XMLObject("String");
		    modelDir.addAttribute("default", application.getAttribute("modelDir"));
		    modelDir.addAttribute("name", "modelDir");
		    modelDir.addAttribute("ownerType", "Application");
		    application.removeAttribute("modelDir");
		    application.addChild(modelDir);

		    XMLObject modelInput = new XMLObject("String");
		    modelInput.addAttribute("default", application.getAttribute("modelInput"));
		    modelInput.addAttribute("name", "modelInput");
		    modelInput.addAttribute("ownerType", "Application");
		    application.removeAttribute("modelInput");
		    application.addChild(modelInput);

		    application.addAttribute("name", "Evaluation");
		    application.addAttribute("altName", "Application");
		    application.addAttribute("base", "EvaluationType");
		    application.addAttribute("ownerType", "ProblemSetup");
		}
		else if (functionType.equals("Library"))
		{
		    XMLObject library = problem.getChild("Library");

		    XMLObject init = new XMLObject("String");
		    init.addAttribute("default", library.getAttribute("Init"));
		    init.addAttribute("name", "Init");
		    init.addAttribute("ownerType", "Library");
		    library.removeAttribute("Init");
		    library.addChild(init);

		    XMLObject fEval = new XMLObject("String");
		    fEval.addAttribute("default", library.getAttribute("FEval"));
		    fEval.addAttribute("name", "FEval");
		    fEval.addAttribute("ownerType", "Library");
		    library.removeAttribute("FEval");
		    library.addChild(fEval);

		    XMLObject libName = new XMLObject("String");
		    libName.addAttribute("default", library.getAttribute("LibName"));
		    libName.addAttribute("name", "LibName");
		    libName.addAttribute("ownerType", "Library");
		    library.removeAttribute("LibName");
		    library.addChild(libName);

		    XMLObject first = new XMLObject("Boolean");
		    first.addAttribute("default", library.getAttribute("First"));
		    first.addAttribute("name", "First");
		    first.addAttribute("ownerType", "Library");
		    library.removeAttribute("First");
		    library.addChild(first);

		    XMLObject second = new XMLObject("Boolean");
		    second.addAttribute("default", library.getAttribute("Second"));
		    second.addAttribute("name", "Second");
		    second.addAttribute("ownerType", "Library");
		    library.removeAttribute("Second");
		    library.addChild(second);

		    library.addAttribute("name", "Evaluation");
		    library.addAttribute("altName", "Library");
		    library.addAttribute("base", "EvaluationType");
		    library.addAttribute("ownerType", "ProblemSetup");
		}

		XMLObject variables = problem.getChild("VariableClass");
		XMLObject lConstraints = problem.getChild("LConstraintClass");
		XMLObject nlConstraints = problem.getChild("NLConstraintClass");

		XMLObject numVars = new XMLObject("Integer");
		numVars.addAttribute("default", variables.getAttribute("numVariables"));
		numVars.addAttribute("name", "numVariables");
		numVars.addAttribute("ownerType", "VariableClass");
		variables.removeAttribute("numVariables");
		variables.addChild(numVars);
		variables.addAttribute("name", "VariablesBounds");
		variables.addAttribute("altName", "VariableClass");
		variables.addAttribute("base", "VariableClass");
		variables.addAttribute("ownerType", "ProblemSetup");

		XMLObject numLCons = new XMLObject("Integer");
		numLCons.addAttribute("default", lConstraints.getAttribute("numLConstraints"));
		numLCons.addAttribute("name", "numLConstraints");
		numLCons.addAttribute("ownerType", "LConstraintClass");
		lConstraints.removeAttribute("numLConstraints");
		lConstraints.addChild(numLCons);
		lConstraints.addAttribute("name", "LinearConstraints");
		lConstraints.addAttribute("altName", "LConstraintClass");
		lConstraints.addAttribute("base", "LConstraintClass");
		lConstraints.addAttribute("ownerType", "ProblemSetup");

		XMLObject numNLCons = new XMLObject("Integer");
		numNLCons.addAttribute("default", nlConstraints.getAttribute("numNLConstraints"));
		numNLCons.addAttribute("name", "numNLConstraints");
		numNLCons.addAttribute("ownerType", "NLConstraintClass");
		nlConstraints.removeAttribute("numNLConstraints");
		nlConstraints.addChild(numNLCons);
		nlConstraints.addAttribute("name", "NonLinearConstraints");
		nlConstraints.addAttribute("altName", "NLConstraintClass");
		nlConstraints.addAttribute("base", "NLConstraintClass");
		nlConstraints.addAttribute("ownerType", "ProblemSetup");

		int numVariables = variables.getChild("Array").numChildren();
		int numLConstraints = lConstraints.getChild("Array").numChildren();
		int numNLConstraints = nlConstraints.getChild("Array").numChildren();

		XMLObject vContents = new XMLObject("Contents");
		XMLObject lcContents = new XMLObject("Contents");
		XMLObject nlcContents = new XMLObject("Contents");
    
		for (int i=0; i<numVariables; i++)
		{
		    XMLObject item = new XMLObject("Item");
		    item.addAttribute("index", String.valueOf(i));

		    XMLObject var_i = variables.getChild("Array").getChild(i);

		    XMLObject theName = new XMLObject("String");
		    theName.addAttribute("default", var_i.getAttribute("theName"));
		    theName.addAttribute("name", "theName");
		    theName.addAttribute("ownerType", "Variables");
		    var_i.removeAttribute("theName");
		    var_i.addChild(theName);

		    XMLObject initVal = new XMLObject("Double");
		    initVal.addAttribute("default", var_i.getAttribute("initVal"));
		    initVal.addAttribute("name", "initVal");
		    initVal.addAttribute("ownerType", "Variables");
		    var_i.removeAttribute("initVal");
		    var_i.addChild(initVal);

		    XMLObject lower = new XMLObject("Double");
		    lower.addAttribute("default", var_i.getAttribute("lower"));
		    lower.addAttribute("name", "lower");
		    lower.addAttribute("ownerType", "Variables");
		    var_i.removeAttribute("lower");
		    var_i.addChild(lower);

		    XMLObject upper = new XMLObject("Double");
		    upper.addAttribute("default", var_i.getAttribute("upper"));
		    upper.addAttribute("name", "upper");
		    upper.addAttribute("ownerType", "Variables");
		    var_i.removeAttribute("upper");
		    var_i.addChild(upper);

		    var_i.addAttribute("name", "Variable");
		    var_i.addAttribute("altName", "Variables");
		    var_i.addAttribute("base", "Variables");
		    var_i.addAttribute("ownerType", "VariableClass");

		    item.addChild(var_i);
		    vContents.addChild(item);
		}

		for (int i=0; i<numVariables; i++)
		{
		    variables.getChild("Array").removeChild("Variables");
		}

		for (int i=0; i<numLConstraints; i++)
		{
		    XMLObject item = new XMLObject("Item");
		    item.addAttribute("index", String.valueOf(i));

		    XMLObject lcons_i = lConstraints.getChild("Array").getChild(i);

		    XMLObject constraintName = new XMLObject("String");
		    constraintName.addAttribute("default", lcons_i.getAttribute("constraintName"));
		    constraintName.addAttribute("name", "constraintName");
		    constraintName.addAttribute("ownerType", "LinearConstraint");
		    lcons_i.removeAttribute("constraintName");
		    lcons_i.addChild(constraintName);

		    XMLObject header = new XMLObject("Header");
		    header.addAttribute("label", "($Coefficient)$Variable");
		    header.addAttribute("name", "Coefficient");

		    XMLObject headerName = new XMLObject("String");
		    headerName.addAttribute("default", "OPT_");
		    headerName.addAttribute("name", "Variable");
		    header.addChild(headerName);

		    XMLObject headerValue = new XMLObject("Double");
		    headerValue.addAttribute("default", "0");
		    headerValue.addAttribute("name", "Coefficient");
		    header.addChild(headerValue);

		    XMLObject entries = new XMLObject("Entries");
		    int numEntries = lcons_i.getChild("Table").numChildren();

		    for (int j=0; j<numEntries; j++)
		    {
			XMLObject entry_j = lcons_i.getChild("Table").getChild(j);
			XMLObject vEntry = new XMLObject("Cell");
			vEntry.addAttribute("value", entry_j.getAttribute("Variable"));
			vEntry.addAttribute("field", "Variable");
			entry_j.removeAttribute("Variable");
			entry_j.addChild(vEntry);

			XMLObject cEntry = new XMLObject("Cell");
			cEntry.addAttribute("value", entry_j.getAttribute("Coefficient"));
			cEntry.addAttribute("field", "Coefficient");
			entry_j.removeAttribute("Coefficient");
			entry_j.addChild(cEntry);

			entries.addChild(entry_j);
		    }

		    for (int j=0; j<numEntries; j++)
		    {
			lcons_i.getChild("Table").removeChild("Entry");
		    }

		    lcons_i.getChild("Table").addChild(header);
		    lcons_i.getChild("Table").addChild(entries);

		    XMLObject operator = new XMLObject("String");
		    operator.addAttribute("default", lcons_i.getAttribute("operator"));
		    operator.addAttribute("name", "operator");
		    operator.addAttribute("ownerType", "LinearConstraint");
		    lcons_i.removeAttribute("operator");
		    lcons_i.addChild(operator);

		    XMLObject rhs = new XMLObject("Double");
		    rhs.addAttribute("default", lcons_i.getAttribute("rhs"));
		    rhs.addAttribute("name", "rhs");
		    rhs.addAttribute("ownerType", "LinearConstraint");
		    lcons_i.removeAttribute("rhs");
		    lcons_i.addChild(rhs);

		    lcons_i.addAttribute("name", "Constraint");
		    lcons_i.addAttribute("altName", "LinearConstraint");
		    lcons_i.addAttribute("base", "LinearConstraint");
		    lcons_i.addAttribute("ownerType", "LConstraintClass");

		    item.addChild(lcons_i);
		    lcContents.addChild(item);
		}

		for (int i=0; i<numLConstraints; i++)
		{
		    lConstraints.getChild("Array").removeChild("LinearConstraint");
		}

		for (int i=0; i<numNLConstraints; i++)
		{
		    XMLObject item = new XMLObject("Item");
		    item.addAttribute("index", String.valueOf(i));

		    XMLObject nlcons_i = nlConstraints.getChild("Array").getChild(i);
		    String nlConstraintType = nlcons_i.getChild(1).getTag();

		    XMLObject constraintName = new XMLObject("String");
		    constraintName.addAttribute("default", nlcons_i.getAttribute("theName"));
		    constraintName.addAttribute("name", "theName");
		    constraintName.addAttribute("ownerType", "NonLinearConstraint");
		    nlcons_i.removeAttribute("theName");
		    nlcons_i.addChild(constraintName);

		    if (nlConstraintType.equals("Application"))
		    {
			XMLObject application = nlcons_i.getChild("Application");

			XMLObject scriptName = new XMLObject("String");
			scriptName.addAttribute("default", application.getAttribute("scriptName"));
			scriptName.addAttribute("name", "scriptName");
			scriptName.addAttribute("ownerType", "Application");
			application.removeAttribute("scriptName");
			application.addChild(scriptName);

			XMLObject modelDir = new XMLObject("String");
			modelDir.addAttribute("default", application.getAttribute("modelDir"));
			modelDir.addAttribute("name", "modelDir");
			modelDir.addAttribute("ownerType", "Application");
			application.removeAttribute("modelDir");
			application.addChild(modelDir);

			XMLObject modelInput = new XMLObject("String");
			modelInput.addAttribute("default", application.getAttribute("modelInput"));
			modelInput.addAttribute("name", "modelInput");
			modelInput.addAttribute("ownerType", "Application");
			application.removeAttribute("modelInput");
			application.addChild(modelInput);

			application.addAttribute("name", "Evaluation");
			application.addAttribute("altName", "Application");
			application.addAttribute("base", "EvaluationType");
			application.addAttribute("ownerType", "NonLinearConstraint");
		    }
		    else if (nlConstraintType.equals("Library"))
		    {
			XMLObject library = nlcons_i.getChild("Library");

			XMLObject init = new XMLObject("String");
			init.addAttribute("default", library.getAttribute("Init"));
			init.addAttribute("name", "Init");
			init.addAttribute("ownerType", "Library");
			library.removeAttribute("Init");
			library.addChild(init);

			XMLObject fEval = new XMLObject("String");
			fEval.addAttribute("default", library.getAttribute("FEval"));
			fEval.addAttribute("name", "FEval");
			fEval.addAttribute("ownerType", "Library");
			library.removeAttribute("FEval");
			library.addChild(fEval);

			XMLObject libName = new XMLObject("String");
			libName.addAttribute("default", library.getAttribute("LibName"));
			libName.addAttribute("name", "LibName");
			libName.addAttribute("ownerType", "Library");
			library.removeAttribute("LibName");
			library.addChild(libName);

			XMLObject first = new XMLObject("Boolean");
			first.addAttribute("default", library.getAttribute("First"));
			first.addAttribute("name", "First");
			first.addAttribute("ownerType", "Library");
			library.removeAttribute("First");
			library.addChild(first);

			XMLObject second = new XMLObject("Boolean");
			second.addAttribute("default", library.getAttribute("Second"));
			second.addAttribute("name", "Second");
			second.addAttribute("ownerType", "Library");
			library.removeAttribute("Second");
			library.addChild(second);

			library.addAttribute("name", "Evaluation");
			library.addAttribute("altName", "Library");
			library.addAttribute("base", "EvaluationType");
			library.addAttribute("ownerType", "NonLinearConstraint");
		    }

		    XMLObject operator = new XMLObject("String");
		    operator.addAttribute("default", nlcons_i.getAttribute("operator"));
		    operator.addAttribute("name", "operator");
		    operator.addAttribute("ownerType", "NonLinearConstraint");
		    nlcons_i.removeAttribute("operator");
		    nlcons_i.addChild(operator);

		    XMLObject rhs = new XMLObject("Double");
		    rhs.addAttribute("default", nlcons_i.getAttribute("rhs"));
		    rhs.addAttribute("name", "rhs");
		    rhs.addAttribute("ownerType", "NonLinearConstraint");
		    nlcons_i.removeAttribute("rhs");
		    nlcons_i.addChild(rhs);

		    nlcons_i.addAttribute("name", "Constraint");
		    nlcons_i.addAttribute("altName", "NonLinearConstraint");
		    nlcons_i.addAttribute("base", "NonLinearConstraint");
		    nlcons_i.addAttribute("ownerType", "NLConstraintClass");

		    item.addChild(nlcons_i);
		    nlcContents.addChild(item);
		}

		for (int i=0; i<numNLConstraints; i++)
		{
		    nlConstraints.getChild("Array").removeChild("NonLinearConstraint");
		}

		variables.getChild("Array").addChild(vContents);
		lConstraints.getChild("Array").addChild(lcContents);
		nlConstraints.getChild("Array").addChild(nlcContents);

		XMLObject basicOptions = solver.getChild("BasicOptions");

		XMLObject outFile = new XMLObject("String");
		outFile.addAttribute("default", basicOptions.getAttribute("outFile"));
		outFile.addAttribute("name", "outFile");
		outFile.addAttribute("ownerType", "BasicOptions");
		basicOptions.removeAttribute("outFile");
		basicOptions.addChild(outFile);

		XMLObject maxIter = new XMLObject("Integer");
		maxIter.addAttribute("default", basicOptions.getAttribute("maxIter"));
		maxIter.addAttribute("name", "maxIter");
		maxIter.addAttribute("ownerType", "BasicOptions");
		basicOptions.removeAttribute("maxIter");
		basicOptions.addChild(maxIter);

		XMLObject maxFeval = new XMLObject("Integer");
		maxFeval.addAttribute("default", basicOptions.getAttribute("maxFeval"));
		maxFeval.addAttribute("name", "maxFeval");
		maxFeval.addAttribute("ownerType", "BasicOptions");
		basicOptions.removeAttribute("maxFeval");
		basicOptions.addChild(maxFeval);

		XMLObject debug = new XMLObject("Boolean");
		debug.addAttribute("default", basicOptions.getAttribute("Debug"));
		debug.addAttribute("name", "Debug");
		debug.addAttribute("ownerType", "BasicOptions");
		basicOptions.removeAttribute("Debug");
		basicOptions.addChild(debug);

		basicOptions.addAttribute("name", "BasicOptions");
		basicOptions.addAttribute("altName", "BasicOptions");
		basicOptions.addAttribute("base", "BasicOptions");
		basicOptions.addAttribute("ownerType", "AlgorithmDef");

		XMLObject advancedOptions = solver.getChild("AdvancedOptions");

		XMLObject fcnTol = new XMLObject("Double");
		fcnTol.addAttribute("default", advancedOptions.getAttribute("fcnTol"));
		fcnTol.addAttribute("name", "fcnTol");
		fcnTol.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("fcnTol");
		advancedOptions.addChild(fcnTol);

		XMLObject gradTol = new XMLObject("Double");
		gradTol.addAttribute("default", advancedOptions.getAttribute("gradTol"));
		gradTol.addAttribute("name", "gradTol");
		gradTol.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("gradTol");
		advancedOptions.addChild(gradTol);

		XMLObject stepTol = new XMLObject("Double");
		stepTol.addAttribute("default", advancedOptions.getAttribute("stepTol"));
		stepTol.addAttribute("name", "stepTol");
		stepTol.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("stepTol");
		advancedOptions.addChild(stepTol);

		XMLObject conTol = new XMLObject("Double");
		conTol.addAttribute("default", advancedOptions.getAttribute("conTol"));
		conTol.addAttribute("name", "conTol");
		conTol.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("conTol");
		advancedOptions.addChild(conTol);

		XMLObject minStep = new XMLObject("Double");
		minStep.addAttribute("default", advancedOptions.getAttribute("minStep"));
		minStep.addAttribute("name", "minStep");
		minStep.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("minStep");
		advancedOptions.addChild(minStep);

		XMLObject maxStep = new XMLObject("Double");
		maxStep.addAttribute("default", advancedOptions.getAttribute("maxStep"));
		maxStep.addAttribute("name", "maxStep");
		maxStep.addAttribute("ownerType", "AdvancedOptions");
		advancedOptions.removeAttribute("maxStep");
		advancedOptions.addChild(maxStep);

		advancedOptions.addAttribute("name", "AdvancedOptions");
		advancedOptions.addAttribute("altName", "AdvancedOptions");
		advancedOptions.addAttribute("base", "AdvancedOptions");
		advancedOptions.addAttribute("ownerType", "AlgorithmDef");

		if (solverType.equals("NIPS"))
		{
		    solver.addAttribute("name", "AlgorithmDefinition");
		    solver.addAttribute("altName", "NIPS");
		    solver.addAttribute("base", "AlgorithmDef");
		    solver.addAttribute("ownerType", "OPT");

		    String meritFunction = solver.getChild(3).getTag();
		    if (meritFunction.equals("argaezTapia"))
		    {
			XMLObject argaezTapia = solver.getChild("argaezTapia");

			XMLObject centParm = new XMLObject("Double");
			centParm.addAttribute("default", argaezTapia.getAttribute("centParm"));
			centParm.addAttribute("name", "centParm");
			centParm.addAttribute("ownerType", "argaezTapia");
			argaezTapia.removeAttribute("centParm");
			argaezTapia.addChild(centParm);

			XMLObject lenBound = new XMLObject("Double");
			lenBound.addAttribute("default", argaezTapia.getAttribute("lenBound"));
			lenBound.addAttribute("name", "lenBound");
			lenBound.addAttribute("ownerType", "argaezTapia");
			argaezTapia.removeAttribute("lenBound");
			argaezTapia.addChild(lenBound);

			XMLObject maxBack = new XMLObject("Integer");
			maxBack.addAttribute("default", argaezTapia.getAttribute("maxBack"));
			maxBack.addAttribute("name", "maxBack");
			maxBack.addAttribute("ownerType", "argaezTapia");
			argaezTapia.removeAttribute("maxBack");
			argaezTapia.addChild(maxBack);

			argaezTapia.addAttribute("name", "NIPSParameters");
			argaezTapia.addAttribute("altName", "argaezTapia");
			argaezTapia.addAttribute("base", "NIPSvar");
			argaezTapia.addAttribute("ownerType", "NIPS");
		    }
		    else if (meritFunction.equals("NormFMu"))
		    {
			XMLObject normFMu = solver.getChild("NormFMu");

			XMLObject centParm = new XMLObject("Double");
			centParm.addAttribute("default", normFMu.getAttribute("centParm"));
			centParm.addAttribute("name", "centParm");
			centParm.addAttribute("ownerType", "NormFMu");
			normFMu.removeAttribute("centParm");
			normFMu.addChild(centParm);

			XMLObject lenBound = new XMLObject("Double");
			lenBound.addAttribute("default", normFMu.getAttribute("lenBound"));
			lenBound.addAttribute("name", "lenBound");
			lenBound.addAttribute("ownerType", "NormFMu");
			normFMu.removeAttribute("lenBound");
			normFMu.addChild(lenBound);

			XMLObject maxBack = new XMLObject("Integer");
			maxBack.addAttribute("default", normFMu.getAttribute("maxBack"));
			maxBack.addAttribute("name", "maxBack");
			maxBack.addAttribute("ownerType", "NormFMu");
			normFMu.removeAttribute("maxBack");
			normFMu.addChild(maxBack);

			normFMu.addAttribute("name", "NIPSParameters");
			normFMu.addAttribute("altName", "NormFMu");
			normFMu.addAttribute("base", "NIPSvar");
			normFMu.addAttribute("ownerType", "NIPS");
		    }
		    else if (meritFunction.equals("VanShanno"))
		    {
			XMLObject vanShanno = solver.getChild("VanShanno");

			XMLObject centParm = new XMLObject("Double");
			centParm.addAttribute("default", vanShanno.getAttribute("centParm"));
			centParm.addAttribute("name", "centParm");
			centParm.addAttribute("ownerType", "VanShanno");
			vanShanno.removeAttribute("centParm");
			vanShanno.addChild(centParm);

			XMLObject lenBound = new XMLObject("Double");
			lenBound.addAttribute("default", vanShanno.getAttribute("lenBound"));
			lenBound.addAttribute("name", "lenBound");
			lenBound.addAttribute("ownerType", "VanShanno");
			vanShanno.removeAttribute("lenBound");
			vanShanno.addChild(lenBound);

			XMLObject maxBack = new XMLObject("Integer");
			maxBack.addAttribute("default", vanShanno.getAttribute("maxBack"));
			maxBack.addAttribute("name", "maxBack");
			maxBack.addAttribute("ownerType", "VanShanno");
			vanShanno.removeAttribute("maxBack");
			vanShanno.addChild(maxBack);

			vanShanno.addAttribute("name", "NIPSParameters");
			vanShanno.addAttribute("altName", "VanShanno");
			vanShanno.addAttribute("base", "NIPSvar");
			vanShanno.addAttribute("ownerType", "NIPS");
		    }
		}
		else if (solverType.equals("PDS"))
		{
		    solver.addAttribute("name", "AlgorithmDefinition");
		    solver.addAttribute("altName", "PDS");
		    solver.addAttribute("base", "AlgorithmDef");
		    solver.addAttribute("ownerType", "OPT");

		    XMLObject pdsVar = solver.getChild("PDSvar");

		    XMLObject simpType = new XMLObject("String");
		    if (pdsVar.getAttribute("simpType").equals("1"))
		    {
			simpType.addAttribute("default", "right angle simplex");
		    }
		    else if (pdsVar.getAttribute("simpType").equals("2"))
		    {
			simpType.addAttribute("default", "regular simplex");
		    }
		    else
		    {
			simpType.addAttribute("default", "scaled right-angle simplex");
		    }
		    simpType.addAttribute("name", "simpType");
		    simpType.addAttribute("ownerType", "PDSvar");
		    pdsVar.removeAttribute("simpType");
		    pdsVar.addChild(simpType);

		    XMLObject searchSize = new XMLObject("Integer");
		    searchSize.addAttribute("default", pdsVar.getAttribute("searchSize"));
		    searchSize.addAttribute("name", "searchSize");
		    searchSize.addAttribute("ownerType", "PDSvar");
		    pdsVar.removeAttribute("searchSize");
		    pdsVar.addChild(searchSize);

		    pdsVar.addAttribute("name", "PDSParameters");
		    pdsVar.addAttribute("altName", "PDSvar");
		    pdsVar.addAttribute("base", "PDSvar");
		    pdsVar.addAttribute("ownerType", "PDS");
		}
		if (solverType.equals("Newton"))
		{
		    solver.addAttribute("name", "AlgorithmDefinition");
		    solver.addAttribute("altName", "Newton");
		    solver.addAttribute("base", "AlgorithmDef");
		    solver.addAttribute("ownerType", "OPT");

		    String searchStrategy = solver.getChild(3).getTag();
		    if (searchStrategy.equals("lineSearch"))
		    {
			XMLObject lineSearch = solver.getChild("lineSearch");

			XMLObject maxBTIter = new XMLObject("Integer");
			maxBTIter.addAttribute("default", lineSearch.getAttribute("maxBTIter"));
			maxBTIter.addAttribute("name", "maxBTIter");
			maxBTIter.addAttribute("ownerType", "lineSearch");
			lineSearch.removeAttribute("maxBTIter");
			lineSearch.addChild(maxBTIter);

			lineSearch.addAttribute("name", "NewtonParameters");
			lineSearch.addAttribute("altName", "lineSearch");
			lineSearch.addAttribute("base", "Newtonvar");
			lineSearch.addAttribute("ownerType", "Newton");
		    }
		    else if (searchStrategy.equals("trustRegion"))
		    {
			XMLObject trustRegion = solver.getChild("trustRegion");

			XMLObject gradMult = new XMLObject("Double");
			gradMult.addAttribute("default", trustRegion.getAttribute("gradMult"));
			gradMult.addAttribute("name", "gradMult");
			gradMult.addAttribute("ownerType", "trustRegion");
			trustRegion.removeAttribute("gradMult");
			trustRegion.addChild(gradMult);

			trustRegion.addAttribute("name", "NewtonParameters");
			trustRegion.addAttribute("altName", "trustRegion");
			trustRegion.addAttribute("base", "Newtonvar");
			trustRegion.addAttribute("ownerType", "Newton");
		    }
		    else if (searchStrategy.equals("trustPDS"))
		    {
			XMLObject trustPDS = solver.getChild("trustPDS");

			XMLObject gradMult = new XMLObject("Double");
			gradMult.addAttribute("default", trustPDS.getAttribute("gradMult"));
			gradMult.addAttribute("name", "gradMult");
			gradMult.addAttribute("ownerType", "trustPDS");
			trustPDS.removeAttribute("gradMult");
			trustPDS.addChild(gradMult);

			XMLObject searchSize = new XMLObject("Integer");
			searchSize.addAttribute("default", trustPDS.getAttribute("searchSize"));
			searchSize.addAttribute("name", "searchSize");
			searchSize.addAttribute("ownerType", "trustPDS");
			trustPDS.removeAttribute("searchSize");
			trustPDS.addChild(searchSize);

			trustPDS.addAttribute("name", "NewtonParameters");
			trustPDS.addAttribute("altName", "trustPDS");
			trustPDS.addAttribute("base", "Newtonvar");
			trustPDS.addAttribute("ownerType", "Newton");
		    }
		}
		else if (solverType.equals("CG"))
		{
		    solver.addAttribute("name", "AlgorithmDefinition");
		    solver.addAttribute("altName", "CG");
		    solver.addAttribute("base", "AlgorithmDef");
		    solver.addAttribute("ownerType", "OPT");
		}

		/* finally, update the editor's contents with the new data */
		Vector failureList = new Vector();
		editor_.updateValue(body, failureList, false);
	 }
	catch(Exception ex)
	 {
		// Display an error message telling the user that 
		// something has gone wrong. (AJR 03-29-02)
		JOptionPane.showMessageDialog(null, 
																	"There was an error updating Maui with the selected file\n"
																	+ "Maui has been updated as much as possible",
																	"Maui Update Error",
																	JOptionPane.ERROR_MESSAGE);
			
	 }
 }


 /**
	* Returns true if valid XML is required.  If Maui XML is
	* being read in, then there is no need for the current XML
	* to be valid, so ReadAction returns false.
	* @return true if valid XML is requried.
	*/
 public boolean requiresValidXML()
 {
	return false;
 }
}		



