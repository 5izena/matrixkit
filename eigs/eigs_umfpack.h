/*
  In this header file, I have defined a simplified function call to
  the ARPACK solver routine for a complex, asymmetric eigenvalue
  problem Av = lv.  The looping procedure and final extraction of
  eigenvalues and vectors is handled automatically.  Most of the
  parameters to the FORTRAN functions are hidden from the user, since
  most of them are determined from user input anyway.
  
  The remaining parameters to the function calls are as follows:
  
    dsaupd(int n, int nev, complex<double> *Evals)
    dsaupd(int n, int nev, complex<double> *Evals, double **Evecs)

    n: the order of the square matrix A
    nev: the number of eigenvalues to be found, starting at the
         bottom.  Note that the highest eigenvalues, or some
	 other choices, can be found.  For now, the choice of
	 the lowest nev eigenvalues is hard-coded.
    Evals: a one-dimensional array of length nev to hold the
           eigenvalues.
    Evecs: a two-dimensional array of size nev by n to hold the
           eigenvectors.  If this argument is not provided, the
	   eigenvectors are not calculated.  Note that the
	   vectors are stored as columns of Evecs, so that the
	   elements of vector i are the values Evecs[j][i].

  The function is overdefined, so that you can call it in either
  fashion and the appropriate version will be run.

  To use these function calls, there must be a function
  defined in the calling program of the form

    av(int n, complex<double> *in, complex<double> *out)

  where the function av finds out = A.in, and n is the order of the
  matrix A.  This function must be defined before the statement that
  includes this header file, as it needs to know about this function.
  It is used in the looping procedure.

  Note that 0 < nev < n-1.

  Scot Shaw
  30 August 1999 */

#ifndef EIGS_H
#define EIGS_H

#ifdef __cplusplus
extern "C" {
#endif

	
/**
*> DLAMCH determines double precision machine parameters.
*> \param[in] CMACH
*> \verbatim
*>          Specifies the value to be returned by DLAMCH:
*>          = 'E' or 'e',   DLAMCH := eps
*>          = 'S' or 's ,   DLAMCH := sfmin
*>          = 'B' or 'b',   DLAMCH := base
*>          = 'P' or 'p',   DLAMCH := eps*base
*>          = 'N' or 'n',   DLAMCH := t
*>          = 'R' or 'r',   DLAMCH := rnd
*>          = 'M' or 'm',   DLAMCH := emin
*>          = 'U' or 'u',   DLAMCH := rmin
*>          = 'L' or 'l',   DLAMCH := emax
*>          = 'O' or 'o',   DLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
*> \endverbatim
*/
extern "C" double dlamch_( char *cmach );

/**
	Largest eigenvalues and eigenvectors of matrix
	Largest magnitude (default)(which[3]="LM").
	Note:
	i_max,j_max:DEL, no useful any more.
*/
int arpack_eigs(int argc, char * argv[], const int n, const int * i_p, const int ilen, 
	const int i_max, const int * j_p, const int jlen, const int j_max, const double * s_p_real, 
	const double * s_p_imag, const int slen, const int nmodes, const double shift, const double tol, 
	const int disp, const bool isreal, const int vlen,const int nblocksize,const int nnumblock,
	double * v_real,double * v_imag, double * d_real,double * d_imag);

#ifdef __cplusplus
}
#endif

#endif