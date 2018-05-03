#ifndef OptGSS_h
#define OptGSS_h

/*----------------------------------------------------------------------
  Copyright (c) 2003
 ----------------------------------------------------------------------*/

 /**
  *
  * @author R.A.Oliva (raoliva@lbl.gov) Lawrence Berkely National Laboratories.
  *
  */

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef WITH_MPI
#include "mpi.h"
#endif

#include "OptDirect.h"
#include "GenSet.h"
#include "Opt_PARAMS.h"

using std::cerr;

namespace OPTPP {

class OptGSS : public OptDirect {

protected:

  NLP0* nlp;	
  ///< Pointer to an NLP0 object

  NLP1* nlp1;	
  ///< Pointer to an NLP1 object

  NEWMAT::ColumnVector X;		
  ///< Best (current) point in search 

  double fX;		
  ///< Value of objective function at X

  NEWMAT::ColumnVector gX;		
  ///< Value of gradient at X (when available)

  double fprev; 
  ///< stores previous best value

  double Delta; 
  ///< Step-length parameter

  double Phi;   
  ///< Steplength expanding parameter (>=1)

  double Theta; 
  ///< Steplength contracting parameter (0<Theta<ThetaMax<1)

  double Delta_tol; 
  ///< Convergence tolerance. Algoritms stops when Theta <= Theta_tol 

  int Iter_max;
  ///< Upper limit on the number of iterations

  bool SearchAll;
  ///< search type flag (true ==> search on all direction) -- now defunct

  bool computeGrad;
  ///< flag to compute gradient after 1st iteration. Used by trustGSS.
  
  GenSetBase* gset; 
  ///< Pointer to an Generating Set object

  NEWMAT::Matrix extras;  
  ///< Extra search directions

  bool extras_srched;
  ///< True when extra directions have been searched
  
  bool printCOPYRIGHT;
  ///< if true copyright header is printed before optimization

  void printHeader();
  ///< Prints header of output file before optimization begins.

  bool printXiter;
  ///< flag for printing (up to 3) components of X during iterations

  bool printGiter;
  ///< flag for printing (up to 3) components of gX during iterations

  int mpi_rank; // also used in serial code

#ifdef WITH_MPI
  int mpi_size; 

  void setpid() {    // Determine MPI size and rank
    //
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);   // C style
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //
    //    mpi_size = Comm::Get_size();          // C++ style
    //    mpi_rank = Comm::Get_rank();
  }
#endif

public:

  void setParams();
  ///< Assign default values to algorithm parameters

  void setParams(OptGSS_params op);
  ///< Assign values to parameters as specified in argument

  //--
  // Constructors
  //--

  OptGSS() : nlp(0), nlp1(0), gset(0)
    { strcpy(method, "Generating Set Search"); setParams(); }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   */
  OptGSS(NLP0* p, GenSetBase* g) : nlp(p), nlp1(0), gset(g)
    { strcpy(method, "Generating Set Search with an NLP0"); setParams(); }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   */
  OptGSS(NLP1* p, GenSetBase* g) : nlp(p), nlp1(p), gset(g)
    { strcpy(method, "Generating Set Search with an NLP1"); setParams(); }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param M matrix with extra search directions
   */
  OptGSS(NLP0* p, GenSetBase* g, NEWMAT::Matrix& M) : 
    nlp(p), nlp1(0), gset(g), extras(M) 
    {
      strcpy(method, "Generating Set Search with an NLP0 & extra directions");
      setParams(); 
    }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param M matrix with extra search directions
   */
  OptGSS(NLP1* p, GenSetBase* g, NEWMAT::Matrix& M) : 
    nlp(p), nlp1(p), Delta(0.0), computeGrad(true), gset(g), extras(M) { 
    strcpy(method, "Generating Set Search with an NLP1 & extra directions"); 
    setParams(); 
  }

  // Constructors with options passed in

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param op GSS algorithmic options
   */
  OptGSS(NLP0* p, GenSetBase* g, OptGSS_params op) : nlp(p), nlp1(0), gset(g)
    { strcpy(method, "Generating Set Search with an NLP0"); setParams(op); }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param op GSS algorithmic options
   */
  OptGSS(NLP1* p, GenSetBase* g, OptGSS_params op) : nlp(p), nlp1(p), gset(g)
    { strcpy(method, "Generating Set Search with an NLP1"); setParams(op); }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param M matrix with extra search directions
   * @param op GSS algorithmic options
   */
  OptGSS(NLP0* p, GenSetBase* g, NEWMAT::Matrix& M, OptGSS_params op) : 
    nlp(p), nlp1(0), gset(g), extras(M) 
    {
      strcpy(method, "Generating Set Search with an NLP0 & extra directions");
      setParams(op); 
    }

  /**
   * @param p nonlinear problem object
   * @param g generating search object
   * @param M matrix with extra search directions
   * @param op GSS algorithmic options
   */
  OptGSS(NLP1* p, GenSetBase* g, NEWMAT::Matrix& M, OptGSS_params op) : 
    nlp(p), nlp1(p), Delta(0.0), computeGrad(true), gset(g), extras(M) { 
    strcpy(method, "Generating Set Search with an NLP1 & extra directions"); 
    setParams(op); 
  }
  
  //--
  // Destructor(s)
  //--
  ~OptGSS(){;}

  //--
  // Attribute access methods
  //--

  /**
   * Let the user set the step size
   */
  void setStepSize(double s) {Delta = s;}

  /**
   * Let the user set the step tolerance
   */
  void setStepTol(double s) {Delta_tol = s;}

  /**
   * Let the user set the step increment parameter
   */
  void setStepInc(double s) { 
    if (s>1) {Phi = s; return;}
    //cerr << "Step increment factor must exceed 1.0\n";
  }

  /**
   * Let the user set the step decrement parameter
   */
  void setStepDec(double s) { 
    if (0<s && s<1) {Theta = s; return;}
   // cerr << "Step decrement factor must be in interval (0,1)\n";
  }

  /**
   * Let the user set the max number of iterations
   */
  void setMaxIter(int s) { Iter_max = s;}

  /**
   * Let the user set the search strategy
   */
  void setFullSearch(bool s) { SearchAll = s;}

  bool extras_searched() { return extras_srched; }

  void setPrintX(bool s) {printXiter = s;}
  void setPrintG(bool s) {printGiter = s;}
  void printCopyright(bool doit) { printCOPYRIGHT = doit; }

  void setComputeGrad(bool s) {computeGrad = s;}


  //--
  // Our internal Optimization methods
  //--

  /**
   * Default Step Increment Method
   */
  void expandStep() { Delta *= Phi;}

  /**
   * Default Step Decrement Method
   */
  void contractStep() { Delta *= Theta; }

  /**
   * Default update of X and FX with new values
   */
  void updateX(double& newfX, NEWMAT::ColumnVector& newX) { 
    fX = newfX;
    X  << newX;
  }

  /**
   * Default Stopping criterion
   */
  bool StopCondition() { return Delta_tol > Delta;}
  int  StepCondition();

  /**
   * Search for improved point; 
   */
  int search();
  //  int search(bool flag=false);
  //  int searchExtras() { return search(true);  }
  //  int searchGenSet() { return search(false); }


  /**
   * Reset parameter values 
   */
  virtual void reset();

  //--
  // Optimization class methods
  //--

  // -- these are virtual in Optimizeclass
  virtual void initOpt();
  virtual void optimize();
  virtual int checkConvg();
  int checkConvg_fcn();
  int checkConvg_grad();

#ifdef DAKOTA_OPTPP
  virtual void printStatus(char * msg) { printStatus(msg, true);}
#else
  virtual void printStatus(char * msg) { printStatus(msg, false);}
  /// Print the status of optimizer at Xc, but not Xc 
#endif

  void printStatus(char * msg, bool printXc);
  /// Print the status of optimizer at Xc, and Xc if 2n arg==true

  //--
  // Define unused virtal methods from base class (OptimizeClass)
  //--
  NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& ) 
    {return NEWMAT::ColumnVector();}

  virtual void acceptStep(int k, int step_type)
    {OptimizeClass::defaultAcceptStep(k, step_type);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  virtual void readOptInput() {;}

  void printIter(int iter, int bp);

#ifdef WITH_MPI
  int getMPIRank() { return mpi_rank; }
  int getMPISize() { return mpi_size; }
  //  int pll_search();
#endif

};

} // namespace OPTPP
#endif
