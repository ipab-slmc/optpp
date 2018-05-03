#ifndef __CBLAS_H__
#define __CBLAS_H__

#ifdef __cplusplus
extern "C" {
#endif

/* BLAS TYPE DEFINITIONS  */

#define daxpy daxpy_
#define dswap dswap_
#define ddot ddot_
#define dscal dscal_
#define dnrm2 dnrm2_
#define dcopy dcopy_

typedef int  Integer;

typedef struct { 
            float real;
            float imag; 
} Complex;

typedef struct {
	    double real;
            double imag; 
} Zomplex;

typedef enum { NoTranspose,
               Transpose,
               ConjugateTranspose } MatrixTranspose;

typedef enum { UpperTriangle,
               LowerTriangle } MatrixTriangle;

typedef enum { UnitTriangular,
               NotUnitTriangular } MatrixUnitTriangular;

typedef enum { LeftSide,
               RightSide } OperationSide;


/********                                                             ********
 ********                    LEVEL 1 BLAS                             ********
 ********                                                             ********/
/*                             Swap two vectors                              *
 *                                x <-> y                                    */

extern void dswap( Integer *n,  double *x, Integer *incx,  double *y, 
                   Integer *incy );

/*                              Scale a vector                               *
 *                               x <- alpha*x                                */

extern void  dscal( Integer *n,  double *alpha,  double *x, Integer *incx );

/*                         Copy one vector to another                        *
 *                                  y <- x                                   */

extern void dcopy( Integer *n,  double *x, Integer *incx,  double *y, 
                   Integer *incy );

/*                 Scale a vector then add to another vector                 *
 *                             y <- alpha*x + y                              */

extern void daxpy( Integer *n,  double *alpha,  double *x, Integer *incx, 
                    double *y, Integer *incy );

/*                         Dot product of two vectors                        * 
 *                                 dot <- xTy                                */

extern double  ddot( Integer *n, double *x, Integer *incx,
                                double *y, Integer *incy );

/*                            2-Norm of a vector                             *
 *                              nrm2 <- ||x||2                               */

extern double  dnrm2( Integer *n,  double *x, Integer *incx );

#ifdef __cplusplus
}
#endif


#endif /* !__CBLAS_H__ */

