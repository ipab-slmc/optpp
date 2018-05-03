
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <iostream>
#include <string>

using namespace std;

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "xercesc/dom/DOM.hpp"
#include "xercesc/parsers/AbstractDOMParser.hpp"
#include "xercesc/util/PlatformUtils.hpp"

#include "Problem.h"
#include "NIPSProblem.h"
#include "NewtonProblem.h"
#include "CGProblem.h"
#include "PDSProblem.h"
#include "NPSOLProblem.h"

namespace OPTPP {

/**
 * This Program will take an XML file which has been output from Maui
 * for OPT++ and start up an optimization run 
 */

int main(int argc, char ** argv)
{
  // First make sure we have a filename.

  if(argc < 2)
    {
//      cerr << "Error, No filename given" << endl;
//      cerr << "Usage: " << argv[0] << " <filename>" << endl;
      return 0;
    }

#ifdef WITH_MPI
  int me;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

  const char* const filename = argv[1];

  // Setup everything to read in the XML file.

  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    cout << "Error during initialization! :\n";
    return 1;
  }

  static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
  DOMImplementation *impl =
    DOMImplementationRegistry::getDOMImplementation(gLS);
  DOMBuilder        *parser =
    ((DOMImplementationLS*)impl)->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);

  DOMNodeList* tmpList;

  DOMElement* solverXML;
  DOMNode* tmpNode;
	
  DOMElement* root;

  // Parse the file, and if there were any errors, quit.
  // get the root element

  try
    {
      root = parser->parseURI(filename)->getDocumentElement();
    }catch(XMLException & ex)
      {
//	cerr << "Error parsing " << filename << "!" << endl;
//	cerr << XMLString::transcode(ex.getMessage()) << endl;
	exit(1);
      }

  OptError error;
  Problem *optimization;

  // Go through all posible problem types and instantiate the
  // appropriate solver class.

  if((tmpList = root->getElementsByTagName(XMLString::transcode("NIPS")))->getLength() > 0)
    {
      tmpNode = tmpList->item(0);
      solverXML = (DOMElement *)tmpNode;
      optimization = new NIPSProblem(solverXML);
    }
  else if((tmpList = root->getElementsByTagName(XMLString::transcode("CG")))->getLength() > 0)
    {
      tmpNode = tmpList->item(0);
      solverXML = (DOMElement *)tmpNode;
      optimization = new CGProblem(solverXML);
    }
  else if((tmpList = root->getElementsByTagName(XMLString::transcode("Newton")))->getLength() > 0)
    {
      tmpNode = tmpList->item(0);
      solverXML = (DOMElement *)tmpNode;
      optimization = new NewtonProblem(solverXML);
    }
  else if((tmpList = root->getElementsByTagName(XMLString::transcode("PDS")))->getLength() > 0)
    {
      tmpNode = tmpList->item(0);
      solverXML = (DOMElement *)tmpNode;
      optimization = new PDSProblem(solverXML);
    }
  else 
    {
      // no solver class was specified

      error.value = -3;
      error.msg = "Cannot find a solver " +
	          string(XMLString::transcode(root->getTagName())) +
	          " in the input XML file";
    }
	
  // Sent the Problem setup information to the class. 

  if((tmpList = root->getElementsByTagName(XMLString::transcode("ProblemSetup")))->getLength() > 0)
    {
      tmpNode = tmpList->item(0);
      DOMElement* appXML = (DOMElement *)tmpNode;
      optimization->SetProblemXML(appXML);
    }
  else 
    {
      // no problem class was specified

      error.value = -3;
      error.msg = "Problem not correctly specified in the input XML file";
    }

  // Optimize the Problem.

  error = optimization->optimize();

  // if there was an error, print out the corresponding message.

  if(error.value < 0)
    {
//      cerr << error.msg << endl;
    }

#ifdef WITH_MPI
  MPI_Finalize();
#endif    

  return 0;
}

} // namespace OPTPP
