#ifndef hockfcns_h
#define hockfcns_h

#include "CompoundConstraint.h"

using namespace OPTPP;

/*
 *
 * Declarations for Hock and Schittkowski's functions
 */

/* Initializer for Problem 1 */
void init_hs1(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs1(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

CompoundConstraint* create_constraint_hs1(int n);

/* Initializer for Problem 2 */
void init_hs2(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

CompoundConstraint* create_constraint_hs2(int n);

/* Initializer for Problem 5 */
void init_hs5(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs5(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

CompoundConstraint* create_constraint_hs5(int n);

/* Initializer for Problem 6 - as appears on Vanderbei's website */
void init_hs6(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs6(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void eqn_hs6(int mode, int n, const NEWMAT::ColumnVector& x, 
        NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs6(int n);

/* Initializer for Problem 7 */
void init_hs7(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs7(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void eqn_hs7(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs7(int n);

/* Initializer for Problem 10 */
void init_hs10(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs10(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs10(int mode, int n, const NEWMAT::ColumnVector& x, 
        NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs10(int n);

/* Initializer for Problem 13 */
void init_hs13(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs13(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs13(int mode, int n, const NEWMAT::ColumnVector& x, 
             NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs13(int n);

/* Initializer for Problem 14 */
void init_hs14(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs14(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs14(int mode, int n, const NEWMAT::ColumnVector& x, 
     NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs14(int n);

/* Initializer for Problem 26 */
void init_hs26(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs26(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void eqn_hs26(int mode, int n, const NEWMAT::ColumnVector& x, 
        NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs26(int n);

/* Initializer for Problem 28 */
void init_hs28(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs28(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

CompoundConstraint* create_constraint_hs28(int n);

/* Initializer for Problem 35 */
void init_hs35(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs35(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs35(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

CompoundConstraint* create_constraint_hs35(int n);

/* Initializer for Problem 65 */
void init_hs65(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs65(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void hs65_2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

void ineq_hs65(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

void ineq_hs65_2(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, 
       OptppArray<NEWMAT::SymmetricMatrix>& H, int& result);

CompoundConstraint* create_constraint_hs65(int n);
CompoundConstraint* create_constraint_hs65_2(int n);

/* Initializer for Problem 77 */
void init_hs77(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hs77(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs77(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs77(int n);

/* Initializer for Problem 78 */
void init_hs78(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with finite difference derivative */
void hs78(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);

void ineq_hs78(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, int& result);

CompoundConstraint* create_constraint_hs78(int n);


/* Hock and Schittkowski with analytic derivative */

void hs78_2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
            NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

void ineq_hs78_2(int mode, int n, const NEWMAT::ColumnVector& x, 
       NEWMAT::ColumnVector& fx, NEWMAT::Matrix& g, 
       OptppArray<NEWMAT::SymmetricMatrix>& H, int& result);

CompoundConstraint* create_constraint_hs78_2(int n);

/* Initializer for Unconstrained Problem */
void init_hsuncon(int n, NEWMAT::ColumnVector& x);

/* Hock and Schittkowski with analytic derivative */
void hsuncon(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
                NEWMAT::ColumnVector& g, int& result);
void hsuncon2(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
              NEWMAT::ColumnVector& g, NEWMAT::SymmetricMatrix& H, int& result);

#endif
