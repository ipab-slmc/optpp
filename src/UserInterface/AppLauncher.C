
#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

using namespace std;

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "xercesc/util/XMLString.hpp"

#include "AppLauncher.h"

using NEWMAT::ColumnVector;

namespace OPTPP {

string AppLauncher::appDir_ = "";

AppLauncher::AppLauncher(DOMElement* appXML, bool createDir)
{
  char command[120];

  appName_ = XMLString::transcode(appXML->getAttribute(XMLString::transcode("scriptName")));
  appInput_ = XMLString::transcode(appXML->getAttribute(XMLString::transcode("modelInput")));
  string tmpDir = XMLString::transcode(appXML->getAttribute(XMLString::transcode("modelDir")));
  if(tmpDir != "")
  {
    appDir_ = tmpDir;

#ifdef WITH_MPI
    int me, error;
    char wkdir[120];

    if (createDir) {
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      sprintf(wkdir, "%s.Proc%d", appDir_.c_str(), me);
      sprintf(command, "/bin/mkdir -p %s", wkdir);
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error making the MPI working"
//	     << "directories" << endl;
	exit(1);
      }
      chdir(wkdir);
      sprintf(command, "/bin/ln -fs %s/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error creating symbolic links"
//	     << "to model-related files." << endl;
	exit(1);
      }
      sprintf(command, "/bin/cp -Rf %s/makecopies/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error copying template"
//	     << "files." << endl;
	exit(1);
      }
    }
#else
    int error;

    if (createDir) {
      chdir(appDir_.c_str());
      sprintf(command, "/bin/cp -Rf %s/makecopies/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error copying template"
//	     << "files." << endl;
	exit(1);
      }
    }
#endif
  }
}

AppLauncher::AppLauncher(DOMElement* appXML, VariableList& variables,
			 bool createDir)
{
  variables_ = & variables;
  char command[120];

  appName_ = XMLString::transcode(appXML->getAttribute(XMLString::transcode("scriptName")));
  appInput_ = XMLString::transcode(appXML->getAttribute(XMLString::transcode("modelInput")));
  string tmpDir = XMLString::transcode(appXML->getAttribute(XMLString::transcode("modelDir")));
  if(tmpDir != "")
  {
    appDir_ = tmpDir;
#ifdef WITH_MPI
    int me, error;
    char wkdir[120];

    if (createDir) {
      MPI_Comm_rank(MPI_COMM_WORLD, &me);
      sprintf(wkdir, "%s.Proc%d", appDir_.c_str(), me);
      sprintf(command, "/bin/mkdir -p %s", wkdir);
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error making the MPI working"
//	     << "directories" << endl;
	exit(1);
      }
      chdir(wkdir);
      sprintf(command, "/bin/ln -fs %s/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error creating symbolic links"
//	     << "to model-related files." << endl;
	exit(1);
      }
      sprintf(command, "/bin/cp -Rf %s/makecopies/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error copying template"
//	     << "files." << endl;
	exit(1);
      }
    }
#else
    int error;

    if (createDir) {
      chdir(appDir_.c_str());
      sprintf(command, "/bin/cp -Rf %s/makecopies/* .", appDir_.c_str());
      error = system(command);
      if(error == -1)
      {
//	cerr << "AppLauncher: There was an error copying template"
//	     << "files." << endl;
	exit(1);
      }
    }
#endif
  }
}


void AppLauncher::init_app(int ndim, ColumnVector& x)
{
  if(ndim != variables_->size())
  {
//    cerr << "Inside AppLauncher::init_app : the number of dimensions,"
//	 << ndim << ", does not match the number of variables given, "
//	 << variables_->size() << endl;

    exit(1);
  }

  ColumnVector initialValues = variables_->GetInitialValues();

  for(int i = 0; i < ndim; i++)
  {
    x(i+1) = initialValues(i+1);
  }
}


void AppLauncher::run_app(int ndim, const ColumnVector& x, double& fx,
			  int& result)
{
  // Dump all of the data to a file

  int error = setupin(ndim, x, appInput_.c_str());
  if (error == -1) {
//    cerr << "There was an error setting up the application input file"
//	 << endl;
    exit(1);
    }

  // system call to run the app/script 
  int pid = fork();

  if(pid == -1) // if there was an error, quit
  {
//    cerr << "There has been an error in fork, stopping." << endl;
    exit(2);
  }
  else if(pid == 0) // if this is the child, call exec
  {
    RunFunctionEvaluation(ndim, x);
  }
  else // if this is the parent, wait for the child to finish
  {
    int status;

    waitpid(pid, &status, 0);

    // wait for completion
    if(!WIFEXITED(status))
    {
//      cerr << "A child process has exited improperly, stopping." << endl;
      exit(2);
    }

    // read in result value and set fx to this value
		
    ifstream fin("fvalue.out");

    fin >> fx;

    fin.close();

    result = NLPFunction;
  }
}

void AppLauncher::run_app_nlncon(int ndim, int nlncons, const ColumnVector& x, ColumnVector& fx, int& result)
{
  // Dump all of the data to a file

  int error = setupin(ndim, x, appInput_.c_str());
  if (error == -1) {
//    cerr << "There was an error setting up the application input file"
//	 << endl;
    exit(1);
    }

  // system call to run the app/script 
  int pid = fork();

  if(pid == -1) // if there was an error, quit
  {
//    cerr << "There has been an error in fork, stopping." << endl;
    exit(2);
  }
  else if(pid == 0) // if this is the child, call exec
  {
    RunFunctionEvaluation(ndim, x);
  }
  else // if this is the parent, wait for the child to finish
  {
    int status;

    waitpid(pid, &status, 0);

    // wait for completion
    if(!WIFEXITED(status))
    {
//      cerr << "A child process has exited improperly, stopping." << endl;
      exit(2);
    }

    // read in result value and set fx to this value
		
    ifstream fin("convalue.out");
    for(int i = 0; i < nlncons; i++)
    {
      fin >> fx(i+1);
    }

    fin.close();

    result = NLPFunction;
  }
}

void AppLauncher::RunFunctionEvaluation(int ndim, const ColumnVector & x)
{
  int error = execl(appName_.c_str(), appName_.c_str(), NULL);
  if(error == -1)
  {
//    cerr << "There was an error running exec" << endl;
    exit(1);
  }
}

int AppLauncher::setupin(int ndim, const ColumnVector& x,
			 const char *fileName)
{
  int index, retcode;
  char line[80], newLine[80], fileTmplt[80];
  const char *pattern;
  string varName;
  FILE *inputFile, *inputTmplt;

  /* Open files */

  sprintf(fileTmplt, "%s.Tmplt", fileName);

  if ((inputTmplt = fopen(fileTmplt, "r")) == NULL ) {
    printf("setupin: No input deck template found\n");
    return(-1);
  }

  if ((inputFile = fopen(fileName, "w")) == NULL ) {
    printf("setupin: I can't open the input file.\n");
    return(-1);
  }

  int count = 0;
  while (fgets(line, 80, inputTmplt) != NULL) {

    bool match = false;
    for (int i=0; i<ndim; i++) {
      varName = variables_->GetVariableName(i);
      if (strstr(line, varName.c_str()) != NULL) {
	match = true;
	pattern = varName.c_str();
	index = i+1;
      }
    }

    /* substitute pattern with value of ith power */

    if (match) {
      retcode = substitute_value(newLine, line, pattern, x(index));

      if (retcode == 0 ) {  /* no errors, go ahead and print out line */
	fprintf(inputFile, "%s", newLine);
	count++;
      }
      else {
	printf("setupin: error in substitute_value()\n");
	fclose(inputTmplt);
	fclose(inputFile);
	return(-1);
      }
    }
    else {  /* normal line, i.e. no matching pattern, echo line to coyote.in */
      fprintf(inputFile, "%s", line);
    }
  }

  /* check to make sure we used up all of the parameters */

  if (count != ndim) {
    printf("setupin: number of parameters found is not equal to the ");
    printf("number of optimization parameters\n");
    printf("setupin: number of parameters found = %d\n", count);
    printf("setupin: number of optimization parameters = %d\n", ndim);
  }

  fclose(inputTmplt);
  fclose(inputFile);

  return(0);
}

int AppLauncher::substitute_value(char *newLine, char *line,
				  const char *pattern, double value)
{
  /*
   * Given a character string "line", a "pattern", and a real value
   * replace the pattern with the value given.
   *
   * On exit, newLine will contain the changed line
   */

  char tmpLine[80];
  int  i, slen, plen;

  plen = strlen( pattern ); 
  slen = strlen( line ); 

  if ( plen <= slen ) {

    for ( i=0; i<slen-plen; i++ )
      if ( strncmp(pattern,&line[i],plen) == 0 ) break;

    strncpy( newLine, line, i );
    newLine[i] = '\0';
    strcat( newLine, " " );
    newLine[i+1] = '\0';
    sprintf( tmpLine, "%24.16e ", value );
    strcat( newLine, tmpLine );
    strcat( newLine, &line[i+plen]);
  }
  else {
//    cerr << "AppLauncher: Substitution pattern longer than line." << endl;
    exit(1);
  }

  return(0);  
}

} // namespace OPTPP
