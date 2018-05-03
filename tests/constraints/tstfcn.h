#ifndef tstfcn_h
#define tstfcn_h

/*
 *
 * Declarations for TWAFER
 */

/* Initializer for TWAFER */

void init_twaf(int n, NEWMAT::ColumnVector& x);

/* TWAFER with no analytic derivative */

void twaf(int n, const NEWMAT::ColumnVector& x, double& fx, int& result);

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

/* Scaled version of Rosenbrock with analytic derivative */
void srosen(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
	    NEWMAT::ColumnVector& g, int& result);

/* Scaled version of Rosenbrock with analytic derivative and Hessian*/

void srosen2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
	     NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

void   init_illum(int n, NEWMAT::ColumnVector& x);

void   illum(int, int, const NEWMAT::ColumnVector&, double&, NEWMAT::ColumnVector&, int&);

void   illum2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
	      NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

void init_erosen (int ndim, NEWMAT::ColumnVector& x);
void erosen(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void erosen1(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, 
	     NEWMAT::ColumnVector&, int& result);

void init_penalty1 (int ndim, NEWMAT::ColumnVector& x);
void penalty1(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_penalty2 (int ndim, NEWMAT::ColumnVector& x);
void penalty2(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_epowell (int ndim, NEWMAT::ColumnVector& x);
void epowell(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_trig (int ndim, NEWMAT::ColumnVector& x);
void trig(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void trig_as_a_constraint(int ndim, const NEWMAT::ColumnVector& x, 
		NEWMAT::ColumnVector& fx, int& result);

void init_gen_brown (int ndim, NEWMAT::ColumnVector& x);
void gen_brown(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_vardim (int ndim, NEWMAT::ColumnVector& x);
void vardim(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_chain_singular (int ndim, NEWMAT::ColumnVector& x);
void chain_singular(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_gen_wood (int ndim, NEWMAT::ColumnVector& x);
void gen_wood(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_chain_wood (int ndim, NEWMAT::ColumnVector& x);
void chain_wood(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_broyden (int ndim, NEWMAT::ColumnVector& x);
void broyden_tridiag(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void broyden1a(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void broyden1b(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void broyden2a(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void broyden2b(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);
void tointbroy(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_gen_cragg_levy (int ndim, NEWMAT::ColumnVector& x);
void gen_cragg_levy(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_toint_trig (int ndim, NEWMAT::ColumnVector& x);
void toint_trig(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_chebyquad (int ndim, NEWMAT::ColumnVector& x);
void chebyquad(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);

void init_nelder (int ndim, NEWMAT::ColumnVector& x);
void nelder(int ndim, const NEWMAT::ColumnVector& x, double& fx, int& result);


#endif
