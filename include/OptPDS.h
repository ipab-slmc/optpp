#ifndef Optpds_h
#define Optpds_h

/*------------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation 
 J.C. Meza, Sandia National Laboratories meza@california.sandia.gov
 ------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "Opt.h"
#include "NLP0.h"
#include "NLP.h"

using std::ofstream;
using std::ostream;

namespace OPTPP {

/**
 * OptDirect is a derived class of OptimizeClass and the base class for direct
 * search methods.  In OPT++, OptGA, a genetic algorithm, and OptPDS, 
 * a parallel direct search method, are examples of direct search methods.  
 */

class OptDirect: public OptimizeClass {
 protected:
 public:
  OptDirect(){}
  OptDirect(int n): OptimizeClass(n){}
  OptDirect(int n, TOLS t): OptimizeClass(n,t){}
  virtual ~OptDirect(){}
  virtual void acceptStep(int, int) = 0;
  virtual void updateModel(int, int, NEWMAT::ColumnVector) = 0;

  virtual int checkConvg() {return 0;}
  virtual void optimize() {}
  virtual void readOptInput() {}
  virtual void reset() {}
};

//----------------------------------------------------------------------
// Parallel Direct Search Method
//----------------------------------------------------------------------

/**
 * OptPDS is an implementation of a derivative-free algorithm for
 * unconstrained optimization.  The search direction is driven solely
 * by the function information.  OptPDS is easy to implement on parallel
 * machines.  A special feature of this approach is the ease with which 
 * algorithms can be generated to take advantage of any number of
 * processors and to addapt to any cost ratio of communication to function
 * evaluation.  
 *
 * For a further description of the parallel direct search methods see
 * J. E. Dennis, Jr. and Virginia Torczon, "Direct Search Methods on
 * Parallel Machines," SIAM J. Optimization, Vol. 1, No. 4,
 * pp. 448--474, November 1991.
 */

class OptPDS: public OptDirect {

protected:
  NLP0* nlp; 			///< Pointer to an NLP0 object
  NEWMAT::Matrix simplex;		///< Convex hull of dimension+1 points
  NEWMAT::ColumnVector vscales;	///< Vector of scale factors to multiply vertices 


  int search_scheme_size;	///< Number of points used in search scheme
  int simplex_type;  		///< Type of simplex chosen by the user
  int reset_param;
  double tr_size; 		///< Trust-region radius
  double simplex_size;

  bool create_scheme_flag, first, trpds;
  char schemefile_name[80];

public:
  OptPDS(){}
  OptPDS(NLP0* p): OptDirect(p->getDim()), nlp(p),  
      simplex(p->getDim(),p->getDim()+1), vscales(p->getDim()),
      search_scheme_size(64), simplex_type(2), 
      reset_param(0), tr_size(0.0), simplex_size(0.0), create_scheme_flag(true),
      trpds(false) { strcpy(method,"PDS");
      strcpy(schemefile_name,"SCHEME");  vscales = 1.0; NEWMAT::ColumnVector x
							  = p->getXc();
      double perturb;
      for (int i=1; i <= p->getDim(); i++) {
	for (int j=1; j <= p->getDim()+1; j++) {
	  simplex(i,j) = x(i);
	}
      }
      for (int i=1; i<= p->getDim(); i++) {
	perturb = x(i)*.01;
	simplex(i,i+1) = x(i) + perturb;
      }
  }
  
  OptPDS(NLP0* p, TOLS t): OptDirect(p->getDim(), t), nlp(p),  
      simplex(p->getDim(),p->getDim()+1), vscales(p->getDim()),
      search_scheme_size(64), simplex_type(2),
      reset_param(0), tr_size(0.0), simplex_size(0.0), create_scheme_flag(true),
      trpds(false) { strcpy(method,"PDS");
      strcpy(schemefile_name,"SCHEME"); vscales = 1.0; NEWMAT::ColumnVector x
							 = p->getXc();
      double perturb;
      for (int i=1; i <= p->getDim(); i++) {
	for (int j=1; j <= p->getDim()+1; j++) {
	  simplex(i,j) = x(i);
	}
      }
      for (int i=1; i<=p->getDim(); i++) {
	perturb = x(i)*.01;
	simplex(i,i+1) = x(i) + perturb;
      }
  }

  virtual ~OptPDS(){}
  virtual NEWMAT::ColumnVector computeSearch(NEWMAT::SymmetricMatrix& ) 
      {return NEWMAT::ColumnVector();}

  virtual void acceptStep(int k, int step_type)
    {defaultAcceptStep(k, step_type);}

  virtual void updateModel(int k, int ndim, NEWMAT::ColumnVector x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  virtual void reset() 
    {int ndim = nlp->getDim(); OptimizeClass::defaultReset(ndim);}

//-------------------------------------------------------------------------
// Accessor Methods
//-------------------------------------------------------------------------

/// Set number of points in the search scheme
  void setSSS(int s)                 {search_scheme_size = s;}

/// Set type of simplex used by the algorithm 
  void setSimplexType(int t)         {simplex_type = t;}

/// Set vertex scale factor
  void setScale(NEWMAT::ColumnVector x)      {vscales = x;};

/// Set simplex used by the algorithm 
  void setSimplex(NEWMAT::Matrix &m)         {simplex = m;};

  void setCreateFlag(bool flag=true) {create_scheme_flag = flag;};

/// Override the default value of filename
  void setSchemeFileName(char *s)    {strcpy(schemefile_name,s);};

/// Set trust-region size 
  void setTRSize(double tr=0)        {tr_size = tr;}

/// Set simplex size 
  void setSimplexSize(double len)    {simplex_size = len;}

/// Set first nonlinear iteration to either true or false
  void setNonIter(bool init=false)   {first = init;}
  void setTRPDS(bool trcon=false)   {trpds = trcon;}

  int getSimplexType()      const {return simplex_type;}
  int getSSS()              const {return search_scheme_size;}
  NEWMAT::ColumnVector& getScale()  {return vscales;};
  bool getCreateFlag()      const {return create_scheme_flag;};
  char *getSchemeFileName() {return schemefile_name;};
  double getTRSize()        {return tr_size;}
  double getSimplexSize()   {return simplex_size;}
  bool getNonIter()         {return first;}
  bool getMethod()          {return trpds;}

// These are defined elsewhere

  void initOpt();		///< Initialize algorithmic parameters
  void optimize();		///< Call the optimization method 
  void readOptInput();		///< Read user-specificed input options
  int checkConvg();	///< Check to see if algorithm satisfies conv. criteria
  void printStatus(char *);	///< Print status of PDS method 
};

int pdsinit(NLP0 *, ostream *, int, int, int *, int *, double,
	    double *, double *, double *, int *, double *, double *,
	    double *, double *, double *, char *, double, int, int,
	    double);

int pdsopt(NLP0 *, ostream *, double *, int *, int, char *, int, int,
	   double, int, int, double, double *, double, int, double *,
	   int *, char *, double, double, double *, int, int, int,
	   double);

int pdswork(NLP0 *, ostream *, ofstream *, int, double, int, int, int *, 
	    double, int, double *, double *, int *, double *, double *,
	    int *, int, double, double *, char *, double, double, int,
	    int, int, double, FILE *);

int pdschk(NLP0 *,int, double *, double *, double, double *, int, double);  

} // namespace OPTPP
#endif
