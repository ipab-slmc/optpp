#ifndef tstfcn_h
#define tstfcn_h

/*
 *
 * Declarations for rosenbrock's functions
 */

/* Initializer for Rosenbrock */

void init_rosen(int n, NEWMAT::ColumnVector& x);

/* Rosenbrock with no analytic derivative */

void rosen0(int n, const NEWMAT::ColumnVector& x, double& fx, int& result);

/* Rosenbrock with analytic derivative */
void rosen(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& g,
	   int& result);

/* Rosenbrock with analytic derivative and Hessian*/

void rosen2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
	    NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

/*
 * Declarations for illum functions
 */
void   init_illum(int n, NEWMAT::ColumnVector& x);

void   illum(int, int, const NEWMAT::ColumnVector&, double&, NEWMAT::ColumnVector&, int&);

void   illum2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
	      NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

/*
 * Declarations for Hock and Schittkowski Problem 65
 */

void init_hs65(int n, NEWMAT::ColumnVector& x);
void hs65(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);
void ineq_hs65(int mode, int n, const NEWMAT::ColumnVector& x, 
		NEWMAT::ColumnVector& cfx, NEWMAT::Matrix& g, int& result);

#endif
