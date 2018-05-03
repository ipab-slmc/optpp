
import Maui.Interface.*;
import XML.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;

/**
 * Brings up a file browser for selection of a filename into
 * which a maui object will be stored.
 *
 */

public class SaveOptAction extends SaveAction
{
    private File lastFile_;

    /**
     * Constructor:  with a label for the button, a hook to the editor, and
     * the file browser that will be used.
     */

    public SaveOptAction()
    {
	super();
    }

    /**
     * Creates a MauiFileSystem object and sends the body to the
     * MauiFileSystem.
     *
     * @param e an ActionEvent
     * @param body an XMLObject that contains the name of the data to
     *             be written.  The MauiFileSystem figures out which
     *             file has been selected
     */

    public void doAction(ActionEvent e, XMLObject body)
    {
	// first do some error checking 

	// check to make sure variables are consistient

	try
	{
	    boolean bounded = false;
	    boolean lconstrained = false;
	    boolean nlconstrained = false;
	    Vector inputErrors = new Vector();
	    XMLObject variables = body.getChild("ProblemSetup").getChild("VariableClass");
	    XMLObject lconstraints = body.getChild("ProblemSetup").getChild("LConstraintClass");
	    XMLObject nlconstraints = body.getChild("ProblemSetup").getChild("NLConstraintClass");

	    editor_.setValue(inputErrors);

	    if (!variables.getAttribute("numVariables").equals(""))
	    {
		int numVariables = variables.getInt("numVariables");
		if (numVariables != variables.getChild("Array").numChildren())
		{
		    inputErrors.add("The specified number of variables does not match the number of variables entered.");
		}
		else
		{
		    XMLObject var_i;
		    for (int i=0; i<numVariables; i++)
		    {
			var_i = variables.getChild("Array").getChild(i);
			if (!var_i.getAttribute("initVal").equals(""))
			{
			    double initVal = var_i.getDouble("initVal");
			    if (!var_i.getAttribute("lower").equals(""))
			    {
				bounded = true;
				if (initVal < var_i.getDouble("lower"))
				{
				    inputErrors.add("Variable " + var_i.getAttribute("theName") + " violates the specified lower bound.");
				}
			    }
			    if (!var_i.getAttribute("upper").equals(""))
			    {
				bounded = true;
				if (initVal > var_i.getDouble("upper"))
				{
				    inputErrors.add("Variable " + var_i.getAttribute("theName") + " violates the specified upper bound.");
				}
			    }
			}
		    }
		}
		if (!lconstraints.getAttribute("numLConstraints").equals(""))
		{
		    int numLConstraints = lconstraints.getInt("numLConstraints");
		    if (numLConstraints != lconstraints.getChild("Array").numChildren())
		    {
			inputErrors.add("The specified number of linear constraints does not match the number of linear constraints entered.");
		    }
		    else if (numLConstraints > 0)
		    {
			lconstrained = true;
			double sum = 0.0;
			double value = 0.0;
			XMLObject lcons_i;
			for (int i=0; i<numLConstraints; i++)
			{
			    lcons_i = lconstraints.getChild("Array").getChild(i);
			    for (int j=0; j<lcons_i.getChild("Table").numChildren(); j++)
			    {
				boolean varFound = false;
				for (int k=0; k<numVariables; k++)
				{
				    if (lcons_i.getChild("Table").getChild(j).getAttribute("Variable").equals(variables.getChild("Array").getChild(k).getAttribute("theName")))
				    {
					value = variables.getChild("Array").getChild(k).getDouble("initVal");
					varFound = true;
				    }
				}
				if (varFound)
				{
				    sum = sum + lcons_i.getChild("Table").getChild(j).getDouble("Coefficient") * value;
				}
				else
				{
				    inputErrors.add("Invalid variable " + lcons_i.getChild("Table").getChild(j).getAttribute("Variable") + " found in linear constraint " + lcons_i.getAttribute("constraintName"));
				}
			    }
			    if (!lcons_i.getAttribute("rhs").equals(""))
			    {
				double rhs = lcons_i.getDouble("rhs");
				if (lcons_i.getAttribute("operator").equals(">="))
				{
				    if (sum < rhs)
				    {
					inputErrors.add("The initial values of the variables violate linear constraint " + lcons_i.getAttribute("constraintName"));
				    }
				}
				else
				{
				    if (sum != rhs)
				    {
					inputErrors.add("The initial values of the variables violate linear constraint " + lcons_i.getAttribute("constraintName"));
				    }
				}
			    }
			}
		    }
		}
	    }

	    if (!nlconstraints.getAttribute("numNLConstraints").equals(""))
	    {
		int numNLConstraints = nlconstraints.getInt("numNLConstraints");
		if (numNLConstraints != nlconstraints.getChild("Array").numChildren())
		{
		    inputErrors.add("The specified number of nonlinear constraints does not match the number of nonlinear constraints entered.");
		}
		else if (numNLConstraints > 0)
		{
		    nlconstrained = true;
		}
	    }

	    if ((bounded) || (lconstrained) || (nlconstrained))
	    {
		if (body.getChild("CG") != null)
		{
		    inputErrors.add("CG for constrained problems is not supported.");
		}
		else if (body.getChild("Newton") != null)
		{
		    if (body.getChild("Newton").getChild("trustRegion") != null)
		    {
			inputErrors.add("Newton with trust-region search for constrained problems is not supported.");
		    }
		    else if ((lconstrained) || (nlconstrained))
		    {
			inputErrors.add("Newton for linearly and/or nonlinearly constrained problems is not supported.");
		    }
		}
	    }

	    if (body.getChild("PDS") != null)
	    {
		XMLObject pdsOptions = body.getChild("PDS").getChild("PDSvar");
		String simplex = pdsOptions.getAttribute("simpType");
		pdsOptions.removeAttribute("simpType");
		if (simplex.equals("right angle simplex"))
		{
		    pdsOptions.addAttribute("simpType", "1");
		}
		else if (simplex.equals("regular simplex"))
		{
		    pdsOptions.addAttribute("simpType", "2");
		}
		else
		{
		    pdsOptions.addAttribute("simpType", "3");
		}
	    }

	    if (inputErrors.size() > 0)
	    {
		String [] options = new String[2];
		StringBuffer errorBuffer = new StringBuffer();

		options[0] = "Save Anyway";
		options[1] = "Cancel";

		for (int i=0; i<inputErrors.size(); i++)
		{
		    errorBuffer.append(i).append(" ").append(inputErrors.elementAt(i)).append("\n");
		}

		int result = 
		    JOptionPane.showOptionDialog(null, 
						 "The following errors were found:\n" + errorBuffer.toString() + "\nWhat do you want to do?\n",
						 "User Input Error!",
						 JOptionPane.YES_NO_CANCEL_OPTION,
						 JOptionPane.WARNING_MESSAGE,
						 null,
						 options,
						 options[1]);

		if(result == JOptionPane.CLOSED_OPTION || result == 1)
		{
		    return;
		}

		if(result == 0)
		{
		    JOptionPane.showMessageDialog(null, "Warning: The file produced by this save operation may be invalid input for OPT++.",
						  "Warning!",
						  JOptionPane.ERROR_MESSAGE);
		    super.doAction(e, body);
		}
	    }
	    else
	    {
		super.doAction(e, body);
	    }
	}
	catch(Exception ex)
	{
	    System.err.println("SaveOptAction - Exception Caught: " + ex.toString());
	    ex.printStackTrace();
	}
    }

    public boolean requiresValidXML()
    {
	return false;
    }

    public boolean requiresVerboseXML()
    {
	return false;
    }
}
