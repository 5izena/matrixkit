/**
  This is an program using the header file
  znaupd.h, an implementation of the ARPACK sparse matrix
  solver routine.  See that header file for more
  documentation on its use.
  
*/

// Begin with some standard include files.

#include <algorithm>
#include <stdio.h>
//#include <fstream.h>
#include <fstream>
//#include <complex.h>
#include <complex>
#include <time.h>

using namespace std;


#include <math.h>
#include <complex>
#include <iostream>


#include "eigs_umfpack.h"

//#include "slu_zdefs.h"
#include "umfpack.h"


//////////////////////////////////////////////////////////////////////////
//				arpack routine
//////////////////////////////////////////////////////////////////////////

extern "C" void znaupd_(int *ido, const char *bmat, int *n, const char *which,
	int *nev, double *tol, complex<double> *resid,
	int *ncv, complex<double> *v, int *ldv,
	int *iparam, int *ipntr, complex<double> *workd,
	complex<double> *workl, int *lworkl,
	double *rwork, int *info);

extern "C" void zneupd_(const int *rvec, const char *howmny, int *select,
	complex<double> *d, complex<double> *z, int *ldz,
	const complex<double> *sigma, complex<double> *workev, const char *bmat,
	int *n, const char *which, int *nev, double *tol,
	complex<double> *resid, int *ncv,
	complex<double> *v, int *ldv, int *iparam,
	int *ipntr, complex<double> *workd,
	complex<double> *workl, int *lworkl,
	double *rwork, int *ierr);

//////////////////////////////////////////////////////////////////////////
//				lapack routine
//////////////////////////////////////////////////////////////////////////

/**
*> ZGETRF computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*/
extern "C" void zgetrf_( int *m, int *n, complex<double> *a, int *lda, int *ipiv, int *info);
/**
*> ZGETRS solves a system of linear equations
*>    A * X = B,  A**T * X = B,  or  A**H * X = B
*> with a general N-by-N matrix A using the LU factorization computed
*/
extern "C" void zgetrs_( char *trans, int *n, int *nrhs, complex<double> *a, int *lda, 
	int *ipiv, complex<double> *b, int *ldb, int *info);
/**
*>
*> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*> overflow.
*/
extern "C" double dlapy2_( double *x, double *y );

//////////////////////////////////////////////////////////////////////////
//				mkl blas routing
//////////////////////////////////////////////////////////////////////////
extern "C" void cblas_zcopy (const int n, const void *x, const int incx, void *y, const int incy);

//////////////////////////////////////////////////////////////////////////
//				blas routine
//////////////////////////////////////////////////////////////////////////
extern "C" void zcopy_(int *n, complex<double> *zx, int *incx, complex<double> *zy, int *incy);
extern "C" void zscal_(int *n, complex<double> *za, complex<double> *zx, int *incx);
extern "C" void zaxpy_(int *n, complex<double> *za, complex<double> *zx, int *incx, complex<double> *zy, int *incy);
extern "C" double dznrm2_(int *n, complex<double> *x, int *incx);

/*
  To use the routines from znaupd.h, we need to define our own matrix
  multiplication routine.  Here I show one method that would exploit
  sparsity by storing only the non-zero matrix elements (though for
  this example we won't actually have a sparse matrix to work with).

  There is a routine for the multiplication, and a global variable T
  that will hold the non-zero matrix elements and their indicies.  The
  global integer tlen will hold the number of non-zero elements.  The
  matrix multiplication function needs to have the form below, and
  must be defined before the file znaupd.h is included.  */

//void av(int n, complex<double> *in, complex<double> *out);

/**
	for umfpack lu
*/
void zcopy_workd_to_b(int n, const complex<double> *workd, double *Bx, double *Bz) 
{
	for (int i = 0; i < n; i++) {
		Bx[i] = workd[i].real();
		Bz[i] = workd[i].imag();
	}
}

void zcopy_x_to_workd(int n, const double *Xx, const double *Xz, complex<double> *workd)
{
	for (int i = 0; i < n; i++) {
		workd[i]._Val[0] = Xx[i];
		workd[i]._Val[1] = Xz[i];
	}
}

/**
qsort, partition for complex 
cout(r) == count(c) == count(val)
*/

void swap_(int *r, int *c, complex<double> *val, int left, int right)
{
	swap(r[left], r[right]);
	swap(c[left], c[right]);
	swap(val[left], val[right]);
}

int partition_(int n, int *r, int *c, 
	complex<double> *val, int left, int right)
{
	int index = left;
	int i = r[index];
	int j = c[index];
	long long idx_col = (long long)j*(long long)n+(long long)i, idx_col_1 = 0;
	int i_1 = 0, j_1 = 0;

	swap_(r, c, val, index, right);
	for (int k = left; k < right; k++) {
		i_1 = r[k];
		j_1 = c[k];
		idx_col_1 = (long long)j_1*(long long)n+(long long)i_1;
		if (idx_col_1 < idx_col) { //ascending order
			swap_(r, c, val, index, k);
			index++;
		}

	}
	swap_(r, c, val, right, index);
	return index;
}

void qsort_(int n, int *r, int *c, 
	complex<double> *val, int left, int right) 
{
	if (left >= right) {
		return;
	}
	int index = partition_(n, r, c, val, left, right);
	qsort_(n, r, c, val, left, index - 1);
	qsort_(n, r, c, val, index + 1, right);
}

void arpack_znaupd(int *Ap, int *Ai, double *Ax, double *Az, double *Control,
					int n, int nev, int ncv, complex<double> sigma, double tol, 
					complex<double> *Evals, complex<double> **Evecs)
{
	//time
	time_t e_time,s_time; 
	s_time = clock();
	
	//for umfpack lu
	void *Symbolic, *Numeric;
	double UMFPack_Info [UMFPACK_INFO]={0.0};
	double rnorm = 0.0;//maxnorm of residual

	int ldv = n;
	int p = min(max(2 * nev, 20), n);
	int lworkl = 2*(3*p*p+5*p);//3*maxncv*maxncv + 5*maxncv;
	//int lworkl = 3*ncv*ncv + 5*ncv;
	ncv = p;
	int iparam[12] = {0};
	int ipntr[15] = {0};

	double *Xx=NULL,*Xz=NULL,*Bx=NULL,*Bz=NULL,*rwork=NULL;
	complex<double> *resid=NULL,*v=NULL,*workd=NULL,
		*workev=NULL,*workl=NULL,*d=NULL;
	Xx = new double[n];if (Xx==NULL) goto __EXIT;
	Xz = new double[n];if (Xz==NULL) goto __EXIT;
	Bx = new double[n];if (Bx==NULL) goto __EXIT;
	Bz = new double[n];if (Bz==NULL) goto __EXIT;
	int *select = new int[ncv];
//	int *ipiv = new int[n];//for zgetrf(lapack routine)
	resid = new complex<double>[n];		if (resid==NULL) goto __EXIT;
	v	  = new complex<double>[2*n*p];	if ( v == NULL ) goto __EXIT;
	workd = new complex<double>[3*n];	if (workd==NULL) goto __EXIT;
	workev= new complex<double>[2*2*p]; if (workev==NULL)goto __EXIT;
	workl = new complex<double>[lworkl];if (workl==NULL) goto __EXIT;
	d	  = new complex<double>[ncv];	if ( d == NULL ) goto __EXIT;
//	complex<double> *ax = new complex<double>[n];
//	complex<double> *mx = new complex<double>[n];

	rwork = new double[p];				if (rwork==NULL) goto __EXIT;
//	double *rd = new double[3*ncv];
	
	memset(Xx, 0, n*sizeof(double));
	memset(Xz, 0, n*sizeof(double));
	memset(Bx, 0, n*sizeof(double));
	memset(Bz, 0, n*sizeof(double));

	memset(select, 0, ncv*sizeof(int));
//	memset(ipiv, 0, n*sizeof(int));
	memset(v, 0, 2*n*p*sizeof(complex<double>));
	memset(workd, 0, 3*n*sizeof(complex<double>));
	memset(workev, 0, 2*2*p*sizeof(complex<double>));
	memset(workl, 0, lworkl*sizeof(complex<double>));
	memset(d, 0, ncv*sizeof(complex<double>));
//	memset(ax, 0, n*sizeof(complex<double>));
//	memset(mx, 0, n*sizeof(complex<double>));
	memset(rwork, 0, p*sizeof(double));
//	memset(rd, 0, 3*ncv*sizeof(double));

	//init resid by random data
	srand((unsigned int) time(0));
	for (int i = 0; i < n; i++) {
		resid[i] = complex<double> ( (float)rand()/(RAND_MAX+1.0), (float)rand()/(RAND_MAX+1.0) );
	}

	int ido = 0;
	char bmat[2] = "G";/* Specifies that the right hand side matrix
					   should be the identity matrix; this makes
					   the problem a standard eigenvalue problem. */

	char which[3] = "LM"; 
						/* Ask for the nev eigenvalues of smallest
						   magnitude.  The possible options are
						   LM: largest magnitude
						   SM: smallest magnitude
						   LA: largest real component LR
						   SA: smallest real compoent SR
						   LI: largest imaginary component
						   SI: smallest imaginary component */

	int mode = 3;

	int info = 0;
	int rvec = 1;

	//check ncv
	if (ncv>n) ncv = n;

	if (n < 3) {
		//!!Notice, be careful in editing, important for Invoker 
		cout  << "[MODE SOLVER ERROR]{ARPACK:N must be at least 3.}\n";
		goto __EXIT;
	}
	if (ncv < 0)
	{
		ncv = nev * 2 + 1;

		if (ncv < 20)
			ncv = 20;

		if (ncv > n - 1)
			ncv = n - 1 ;
	} else {

		if (ncv < nev+1 || ncv > n) {
			//!!Notice, be careful in editing, important for Invoker 
			cout << "[MODE SOLVER ERROR]{ARPACK:Invalid number of ncv (must be NEV+1<=NCV and NCV<=N).}\n";
			goto __EXIT;
		}
	}
	if (nev <= 0 || nev >= n - 1) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Invalid number of eigenvalues to extract (must be 0 < k < n-1).}\n";
		goto __EXIT;
	}
	if (ncv <= nev || ncv > n) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Opts.p must be greater than k and less than n.}\n";
		goto __EXIT;
	}

#if _DEBUG
	cout << "n	= " << n << endl;
	cout << "nev	= " << nev << endl;
	cout << "ncv	= " << ncv << endl;
#endif

	//////////////////////////////////////////////////////////////////////////
	int maxitr = max(300, (int)ceil((double)2*n/max(p,1)));
	int ishfts = 1;
	iparam[0] = ishfts;
	iparam[1] = 0;
	iparam[2] = maxitr ;
	iparam[3] = 0;		// NB blocksize in recurrence
	iparam[4] = 0;
	iparam[5] = 0;
	iparam[6] = mode;
	iparam[7] = 0;
	iparam[8] = 0;
	iparam[9] = 0;
	iparam[10] = 0;

	/**
	*
	*     %----------------------------------------------------%
	*     | Factor C (C = A - SIGMA*I) in real arithmetic using	|
	*	  | umfpack subroutine umfpack_zi_symbolic and			|
	*	  | umfpack_zi_numeric.			  					    |
	*     %----------------------------------------------------%
	*	
	*/

	/**
	*	begin umfpack lu decompose
	*/
	e_time = clock();
	//cout<< "ARPACK:eigs, before lu decompose, time: "<<s_time<<" , end:"<<e_time<<" , end-time="<<e_time - s_time<<"ms sec="<<(e_time - s_time)/1000<<"s\n";
	s_time = e_time;

	/* ---------------------------------------------------------------------- */
	/* symbolic factorization */
	/* ---------------------------------------------------------------------- */
	int status = 0;
	status = umfpack_zi_symbolic (n, n, Ap, Ai, Ax, Az, &Symbolic,
		Control, UMFPack_Info) ;
	if (status < 0)
	{
		umfpack_zi_report_info (Control, UMFPack_Info) ;
		umfpack_zi_report_status (Control, status) ;
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Symbolic factorization failed.}" << endl;
		goto __EXIT;
	}

	/* print the symbolic factorization */

	//printf ("\nSymbolic factorization of A: ") ;
	//(void) umfpack_zi_report_symbolic (Symbolic, Control) ;

	// Factor the matrix.
	/* ---------------------------------------------------------------------- */
	/* numeric factorization, matrix A = A-sigma*I */
	/* ---------------------------------------------------------------------- */
	status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
		Control, UMFPack_Info) ;
	if (status < 0)
	{
		umfpack_zi_report_info (Control, UMFPack_Info) ;
		umfpack_zi_report_status (Control, status) ;
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Numeric factorization failed.}" << endl ;
		goto __EXIT;
	}

#if _DEBUG
	printf ("\nFreeing symbolic object:\n") ;
#endif

	umfpack_zi_free_symbolic (&Symbolic) ;

	e_time = clock();
	//cout<<"ARPACK:eigs, lu decompose, time: "<<s_time<<" , end:"<<e_time<<" , end-time="<<e_time - s_time <<"ms sec="<<(e_time - s_time)/1000<<"s\n";
	s_time = e_time;

	/* print the numeric factorization */
	//printf ("\nNumeric factorization of A: ") ;
	//(void) umfpack_zi_report_numeric (Numeric, Control) ;

	//////////////////////////////////////////////////////////////////////////
	//get default tolerance if not set
	if (abs(tol) <= 0) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Invalid tol.}" << endl;
		char cmach[2] = "E";
		tol = dlamch_(cmach);
	}
	int nrhs = 1;
	int ierr;
	info = 1;//init resid by random data, if info=1

	int tick_count = 0, nzupd_time = 0;
	time_t nzupd_time_s = 0, nzupd_time_e = 0;
	nzupd_time_s = clock();

	do {
		znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, 
			&ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, rwork, &info);
		nzupd_time_e = clock();
		nzupd_time += nzupd_time_e - nzupd_time_s;
		nzupd_time_s = nzupd_time_e;
		if (ido == -1) {
			/* ---------------------------------------------------------------------- */
			/* solve Ax=b */
			/* ---------------------------------------------------------------------- */
			//copy from workd[ipntr[0] - 1].real to Bx
			//copy from workd[ipntr[0] - 1].imag to Bz
			zcopy_workd_to_b(n, &workd[ipntr[0] - 1], Bx, Bz);
			status = umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz,
				Numeric, Control, UMFPack_Info);
			//copy from Xx to workd[ipntr[1] - 1].real 
			//copy from Xz to workd[ipntr[1] - 1].imag
			zcopy_x_to_workd(n, Xx, Xz, &workd[ipntr[1] - 1]);

			//  umfpack_zi_report_status (Control, status) ;
			if (status < 0)	{
				//!!Notice, be careful in editing, important for Invoker 
				cout  << "[MODE SOLVER ERROR]{ARPACK:Solve failed, ido=-1.}" << endl;
				goto __EXIT;
			}
			//printf ("\nx (solution of Ax=b): ") ;
			//(void) umfpack_zi_report_vector (n, x, Control) ;
			//rnorm = resid (n, r, b, &workd[ipntr[1] - 1], 0, Ap, Ai, Ax) ;
			//printf ("maxnorm of residual: %g\n\n", rnorm) ;

		} else if (ido == 1) {
			/**
			*
			*           %-------------------------------------------%
			*           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   |
			*           | The user should supply his/her own linear |
			*           | system solver here that takes             |
			*           | workd(ipntr(1)) as the input, and returns |
			*           | the result to workd(ipntr(2)).            |
			*           %-------------------------------------------%
			*
			*/

			zcopy_workd_to_b(n, &workd[ipntr[2] - 1], Bx, Bz);
			status = umfpack_zi_solve (UMFPACK_A, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz,
				Numeric, Control, UMFPack_Info) ;
			zcopy_x_to_workd(n, Xx, Xz, &workd[ipntr[1] - 1]);
			//umfpack_di_report_info (Control, Info) ;
			//umfpack_zi_report_status (Control, status) ;
			if (status < 0)	{
				cout << "[MODE SOLVER ERROR]{ARPACK:Solve failed, ido = 1.}" << endl ;
				goto __EXIT;
			}
		} else if (ido == 2) {
			/**
			*
			*           %---------------------------------------------%
			*           |          Perform  y <--- M*x                |
			*           | Need matrix vector multiplication routine   |
			*           | here that takes workd(ipntr(1)) as input    |
			*           | and returns the result to workd(ipntr(2)).  |
			*           %---------------------------------------------%
			*
			*/
			//mv(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
			int incx = 1, incy = 1;
//			zcopy_( &n, &workd[ipntr[0] - 1], &incx, &workd[ipntr[1] - 1], &incy);
			cblas_zcopy( n, &workd[ipntr[0] - 1], incx, &workd[ipntr[1] - 1], incy );
//			memcpy(&workd[ipntr[1] - 1], (complex<double> *)rhs, n*sizeof(complex<double>));
		} else if (ido == 99) {
			//cout << " ARPACK:eigs has converged.\n";
		} else  {
			//!!Notice, be careful in editing, important for Invoker 
			cout << " [MODE SOLVER ERROR]{ARPACK:Unknown Ido.}\n";
		}
		tick_count++;
		/**
		*           %-----------------------------------------%
		*           | L O O P   B A C K to call DNAUPD again. |
		*           %-----------------------------------------%		
		*/
		continue;
	} while ((ido==1)||(ido==-1)||(ido==2));

	/* We don't need the Numeric object any more, so free it. */
	umfpack_zi_free_numeric (&Numeric) ;
#ifdef _DEBUG
	cout << " _naupd ido cout: " << tick_count << endl;
	cout << " _naupd time : " << nzupd_time << endl;
#endif	

	e_time = clock();
	//cout<<"ARPACK:eigs, naupd_, time: "<<s_time<<" , end:"<<e_time<<" , end-time="<<e_time - s_time<<"ms sec="<<(e_time - s_time)/1000<<"s\n";
	s_time = e_time;

	if (info<0) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Error with znaupd, info = " << info << ".}\n";
		goto __EXIT;
	} else {
		/**
		*
		*        %-------------------------------------------%
		*        | No fatal errors occurred.                 |
		*        | Post-Process using ZNEUPD .                |
		*        |                                           |
		*        | Computed eigenvalues may be extracted.    |
		*        |                                           |
		*        | Eigenvectors may also be computed now if  |
		*        | desired.  (indicated by rvec = .true.)    |
		*        %-------------------------------------------%
		*
		*/
		rvec = 1;
		zneupd_(&rvec, "A", select, d, v, &ldv, &sigma, workev,
			bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);
		/**
		*
		*        %----------------------------------------------%
		*        | Eigenvalues are returned in the one          |
		*        | dimensional array D.  The corresponding      |
		*        | eigenvectors are returned in the first NCONV |
		*        | (=IPARAM(5)) columns of the two dimensional  |
		*        | array V if requested.  Otherwise, an         |
		*        | orthogonal basis for the invariant subspace  |
		*        | corresponding to the eigenvalues in D is     |
		*        | returned in V.                               |
		*        %----------------------------------------------%
		*
		*/
		if (ierr!=0) {
			//!!Notice, be careful in editing, important for Invoker 
			cout << "[MODE SOLVER ERROR]{ARPACK:Error with zneupd, info = " << ierr << ".}\n";
			goto __EXIT;
		} else if (info==1) {
			cout << "Maximum number of iterations reached.\n\n";
		} else if (info==3) {
			cout << "No shifts could be applied during implicit\n";
			cout << "Arnoldi update, try increasing NCV.\n\n";
		}

		e_time = clock();
		//cout<<"ARPACK:eigs, neupd_, time: "<<s_time<<" , end:"<<e_time<<" , end-time="<<e_time - s_time<<"ms sec="<<(e_time - s_time)/1000<<"s\n";
		s_time = e_time;


		/**
		*              %---------------------------%
		*              | Compute the residual norm |
		*              |                           |
		*              |   ||  A*x - lambda*x ||   |
		*              |                           |
		*              | for the NCONV accurately  |
		*              | computed eigenvalues and  |
		*              | eigenvectors.  (iparam(5) |
		*              | indicates how many are    |
		*              | accurate to the requested |
		*              | tolerance)                |
		*              %---------------------------%
		*/
		//ierr == 0
//do not compute the residual norm
/*
		{
			int nconv = iparam[4];
			for (int j = 1; j <= nconv; j++) {
				//int t = j << 8;
				int t = (j - 1)*maxncv;
				av(n, &v[t], ax);
				//mv(n, &v[t - 256], mx);
				//same as mv
				int incx = 1, incy = 1;
				zcopy_( &n, &v[t], &incx, mx, &incy);

				complex<double> d_1 = -d[j-1];
				//zaxpy_(&n, &d_1, &v[t - 256], &incx, ax, &incy);
				zaxpy_(&n, &d_1, mx, &incx, ax, &incy);

				rd[j - 1] = d[j-1].real();
				rd[j + ncv - 1] = d[j - 1].imag();
				rd[j + 2*ncv - 1] = dznrm2_(&n, ax, &incx);

				//rd[j + 2*ncv - 1] /= sqrtl(rd[j - 1]*rd[j - 1]+rd[j + ncv - 1]*rd[j + ncv - 1]);
				//same as sqrtl
				rd[j + 2*ncv - 1] /= dlapy2_(&rd[j - 1], &rd[j + ncv - 1]);
			}
		}
*/
		int i, j;
		for (i=0; i<nev; i++) { Evals[i] = d[i]; /*cout << Evals[i] << endl;*/};
		for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[j][i] = v[i*n+j];


		//sort by 'LM' Largest magnitude
		
		complex<double> temp;
		for (i=0; i<nev; i++) {
			for (j=i; j<nev; j++) {
				if (Evals[j].real() > Evals[i].real()) {
					temp = Evals[j];
					Evals[j] = Evals[i];
					Evals[i] = temp;
					for (int k=0; k<n; k++) {
						temp = Evecs[k][i];
						Evecs[k][i] = Evecs[k][j];
						Evecs[k][j] = temp;
					}
				}
			}
		}
		

	}
	

__EXIT:
	delete []select;
	delete []resid;
	delete []v;
	delete []workd;
	delete []workev;
	delete []d;
//	delete []ax;
//	delete []mx;
	delete []rwork;
	delete []workl;
//	delete []rd;
	delete []Xx;
	delete []Xz;
	delete []Bx;
	delete []Bz;


}

void init_val(double *val)
{
	//val
	//////////////////////////////////////////////////////////////////////////
	val[0] = 193.681904581396;
	val[1] = 194.391406843839;
	val[2] = 194.760774213429;
	val[3] = 193.742062429795;
	val[4] = 194.451564692238;
	val[5] = 194.820932061827;
	val[6] = 192.560008307518;
	val[7] = 193.269510569961;
	val[8] = 193.638877939551;
	val[9] = 191.579015876904;
	val[10] = 192.288518139347;
	val[11] = 192.657885508937;
	val[12] = 180.413577917807;
	val[13] = 181.123080180251;
	val[14] = 181.49244754984;
	val[15] = 190.113383375563;
	val[16] = 14.4538681223551;
	val[17] = 14.6133558441949;
	val[18] = 13.5344862121621;
	val[19] = 14.2439884746055;
	val[20] = 14.6133558441949;
	val[21] = 0.897737556561086 ;
	val[22] = -0.110344827586207;
	val[23] = 0.897737556561086;
	val[24] = -0.110344827586207 ;
	val[25] = 0.897737556561086;
	val[26] = -0.110344827586207 ;
	val[27] = 0.897737556561086;
	val[28] = -0.110344827586207;
	val[29] = 0.897737556561086 ;
	val[30] = -0.110344827586207 ;
	val[31] = 0.197655811756285;
	val[32] = -0.0828083018827952;
	val[33] = 0.897737556561086;
	val[34] = -0.110344827586207;
	val[35] = 0.298580121703854 ;
	val[36] = -0.105109489051095;
	val[37] = 0.298580121703854;
	val[38] = -0.105109489051095;
	val[39] = 0.298580121703854;
	val[40] = -0.105109489051095 ;
	val[41] = 0.298580121703854 ;
	val[42] = -0.105109489051095 ;
	val[43] = 0.298580121703854 ;
	val[44] = -0.105109489051095 ;
	val[45] = 0.7237896980219 ;
	val[46] = -0.105109489051095;
	val[47] = 0.298580121703854 ;
	val[48] = -0.105109489051095;
	val[49] = -0.0727221647728292;
	val[50] = -0.0727221647728292 ;
	val[51] = -0.0727221647728292;
	val[52] = -0.0830382106244175 ;
	val[53] = -0.0830382106244175;
	val[54] = -0.0830382106244175;
	val[55] = 0.746341463414634 ;
	val[56] = 0.746341463414634 ;
	val[57] = 0.746341463414634;
	val[58] = 0.996923076923077;
	val[59] = 0.996923076923077;
	val[60] = 0.996923076923077 ;
	val[61] = 1 ;
	val[62] = 1 ; 
	val[63] = 1  ;  
	val[64] = 1   ; 
	val[65] = 1   ; 
	val[66] = 1   ;   
	val[67] = -0.104398586168731 ;
	val[68] = -0.104398586168731 ;
	val[69] = -0.104398586168731 ;
	val[70] = 0.248275862068966  ;
	val[71] = 0.248275862068966 ;
	val[72] = 0.248275862068966  ;
	val[73] = 0.978686679174484 ;
	val[74] = 0.978686679174484 ;
	val[75] = 0.978686679174484  ;
	val[76] = 1   ;   
	val[77] = 1  ; 
	val[78] = 1  ;  
	val[79] = 1  ;
	val[80] = 1  ; 
	val[81] = 1  ;
	val[82] = 1  ; 
	val[83] = 1  ;   
	val[84] = 1 ; 
	
	//////////////////////////////////////////////////////////////////////////
}

void main()
{
  int n;
  double shift = 188.689270916658;
  n = 21; // The order of the matrix
  int nmodes = 1;
  int slen = 85;
  int vlen = 21;
  double tol = 1.0e-10;
  /*
    Here we generate the matrix T for the multiplication.  It helps if
    we know the number of non-zero matrix elements before we find
    them, so that we can save on storage space.  If not, generate
    enough storage space for T to cover an upper limit.  */

  //row
  int  i_p[85] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,  //iall
	  2  ,   3  ,   5  ,   6   ,  8   ,  9  ,  11  ,  12 ,   14 ,   15 ,   17, 18  ,  20  ,  21 ,//ie
	  1   ,  2   ,  4 ,    5  ,   7 ,    8 ,   10  ,  11 ,   13  ,  14  ,  16,17 ,   19  ,  20, //iw
	  4  ,   5    , 6  ,   7   ,  8   ,  9   , 10  ,  11 ,   12 ,   13,    14, 15 ,   16 ,   17,    18  ,  19  ,  20  ,  21, //in
	  1  ,   2  ,   3  ,   4   ,  5  ,   6   ,  7  ,   8,     9 ,   10 ,   11,12 ,   13 ,   14,    15 ,   16 ,   17  ,  18 //is
  };
  //column
  int  j_p[85] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,  //iall
	  1   ,  2   ,  4 ,    5  ,   7 ,    8 ,   10  ,  11 ,   13  ,  14  ,  16,17 ,   19  ,  20, //iw
	  2  ,   3  ,   5  ,   6   ,  8   ,  9  ,  11  ,  12 ,   14 ,   15 ,   17, 18  ,  20  ,  21 ,//ie
	  1  ,   2  ,   3  ,   4   ,  5  ,   6   ,  7  ,   8,     9 ,   10 ,   11,12 ,   13 ,   14,    15 ,   16 ,   17  ,  18, //is
	  4  ,   5    , 6  ,   7   ,  8   ,  9   , 10  ,  11 ,   12 ,   13,    14, 15 ,   16 ,   17,    18  ,  19  ,  20  ,  21 //in
  };

  double s_p_real[85] = {0};
  double s_p_imag[85] = {0};
   
  init_val(s_p_real);

  //time
  time_t end,time; 
  time = clock();

  cout << "+++++++++++++++++++++ARPACK:eigs--->begin+++++++++++++++++++" << endl;
  int nev, ncv, status;
  std::complex<double> sigma = std::complex<double>(shift, 0.0);


  double Control [UMFPACK_CONTROL];
  int nz = slen, nz1 = 0;
  int *Ap = NULL, *Ai = NULL;
  double *Ax = NULL, *Az = NULL;

  int *Arow = new int[nz];
  int *Acol = new int[nz];
  double *Aval = new double[nz];
  double *Avalz = new double[nz];

  for (int i = 0; i < nz; i++) {
	  Arow[i] = i_p[i] - 1;
	  Acol[i] = j_p[i] - 1;
	  if (Arow[i] == Acol[i]) {
		  Aval[i] = s_p_real[i] - sigma.real();
		  Avalz[i] = s_p_imag[i] - sigma.imag();
	  } else {
		  Aval[i] = s_p_real[i];
		  Avalz[i] = s_p_imag[i];
	  }
  }

  /* get the default control parameters */
  umfpack_zi_defaults (Control) ;

  /* change the default print level for this demo */
  /* (otherwise, nothing will print) */
  Control [UMFPACK_PRL] = 6 ;

  /* print the license agreement */
  //umfpack_zi_report_status (Control, UMFPACK_OK) ;
  Control [UMFPACK_PRL] = 5 ;

  /* print the triplet form of the matrix */
  printf ("\nA: ") ;
  (void) umfpack_zi_report_triplet (n, n, nz, Arow, Acol, Aval, Avalz, Control) ;

  /* convert to column form */
  nz1 = max (nz,1) ;	/* ensure arrays are not of size zero. */
  Ap = new int[(n+1) * sizeof (int)];
  Ai = new int [nz1 * sizeof (int)];
  Ax = new double[nz1 * sizeof (double)];
  Az = new double[nz1 * sizeof (double)];
  if (!Ap || !Ai || !Ax || !Az) {
	  cout << "out of memory" << endl;
	  goto __Exit;
  }

  status = umfpack_zi_triplet_to_col (n, n, nz, Arow, Acol, Aval, Avalz,
	  Ap, Ai, Ax, Az, (int *) NULL) ;

  if (status < 0) {
	  umfpack_zi_report_status (Control, status) ;
	  cout << "umfpack_zi_triplet_to_col failed" << endl;
  }

  /* print the column-form of A */
  printf ("\nA: ") ;
  (void) umfpack_zi_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

  //////////////////////////////////////////////////////////////////////////

  nev = nmodes; // The number of values to calculate
  ncv = vlen;

  std::complex<double> *Evals = new std::complex<double>[nev];
  std::complex<double> **Evecs = new std::complex<double>*[n];
  for (int i=0; i<n; i++) Evecs[i] = new std::complex<double>[nev];

  end = clock();
  cout<<"ARPACK:eigs, init, time: "<<time<<" , end:"<<end<<" , end-time="<<end - time<<"ms sec="<<(end - time)/1000<<"s\n\n";
  time = end;

  cout << "\n+++++++++++++++++++++arpack_znaupd---->begin+++++++++++++++++++" << endl;

  arpack_znaupd(Ap, Ai, Ax, Az, Control, n, nev, ncv, sigma, tol, Evals, Evecs);

  cout << "\n+++++++++++++++++++++arpack_znaupd---->end++++++++++++++++++++" << endl;

  end = clock();
  cout<<"ARPACK:eigs, coveraged time: "<<time<<" , end:"<<end<<" , end-time="<<end - time<<"ms sec="<<(end - time)/1000<<"s\n\n";


  printf("shift: %e  \n", shift);

  for (int i = 0; i < nev; i++) {
	  printf("eigenvalue: [%e,%e]\n", Evals[i]._Val[0], Evals[i]._Val[1]);
	  cout << "eigenvect: "<< endl;
	  for (int j = 0; j < ncv; j++) { 
		  if (j < 10)
			  cout << "[" << Evecs[j][i]._Val[0] << "," << Evecs[j][i]._Val[1] <<"]";
	  }
  }
  cout << endl;

__Exit:
  delete []Ap;
  delete []Ai;
  delete []Ax;
  delete []Az;
  delete []Arow;
  delete []Acol;
  delete []Aval;
  delete []Avalz;

  for (int i=0; i<n; i++) { delete []Evecs[i]; }
  delete []Evecs;

  cout << "+++++++++++++++++++++ARPACK:eigs---->end++++++++++++++++++++" << endl;
  
}

//no implementation
void checkInputs()
{
	//isrealprob
	//issysA = false
	//SquareMatrix=true
	//nmodes
	//sigma
	//whitch='LM'
	//eps= 2.220446049250313e-16, default by the LAPACK auxiliary subroutine dlamch
	//init resid
	//stype = "G"
}

//the entry function
int arpack_eigs(int argc, char * argv[], const int n, 
	const int * i_p, const int ilen, const int i_max, const int * j_p, 
	const int jlen, const int j_max, const double * s_p_real, 
	const double * s_p_imag, const int slen, const int nmodes, const double shift, 
	const double tol, const int disp, const bool isreal, const int vlen,
	const int nblocksize,const int nnumblock, double * v_real,double * v_imag,
	double * d_real,double * d_imag)
{
	//time
	time_t end,time; 
	time = clock();

	//test use [21,21] matrix
	//main_test(); return 0;

	//cout << "+++++++++++++++++++++ARPACK:eigs--->begin+++++++++++++++++++" << endl;
	int nev, ncv, status;

	std::complex<double> sigma = std::complex<double>(shift, 0.0);

	checkInputs();

	double Control [UMFPACK_CONTROL];
	int nz = slen, nz1 = 0;
	std::complex<double> *Evals = NULL;
	std::complex<double> **Evecs = NULL;
	int *Ap = NULL, *Ai = NULL;
	double *Ax = NULL, *Az = NULL;

	int *Arow = new int[nz];
	int *Acol = new int[nz];
	double *Aval = new double[nz];
	double *Avalz = new double[nz];
	if (!Arow || !Acol || !Aval || !Avalz) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Out of memory.}" << endl;
		goto __Exit;
	}

	for (int i = 0; i < nz; i++) {
		//Arow[i] = i_p[i] - 1;
		//Acol[i] = j_p[i] - 1;
		Arow[i] = i_p[i];
		Acol[i] = j_p[i];

		if (Arow[i] == Acol[i]) {
			Aval[i] = s_p_real[i] - sigma.real();
			Avalz[i] = s_p_imag[i] - sigma.imag();
		} else {
			Aval[i] = s_p_real[i];
			Avalz[i] = s_p_imag[i];
		}
#ifdef _DEBUG__
		cout <<"("<<Arow[i]<<","<<Acol[i]<<"){"<<Aval[i]<<","<<Avalz[i]<<"}"<<"	";
#endif
	}
#ifdef _DEBUG__
	cout << endl;
#endif
	/* get the default control parameters */
	umfpack_zi_defaults (Control) ;

	/* change the default print level for this demo */
	/* (otherwise, nothing will print) */
	Control [UMFPACK_PRL] = 6 ;

	/* print the license agreement */
	//umfpack_zi_report_status (Control, UMFPACK_OK) ;
	Control [UMFPACK_PRL] = 5 ;
	
	/* print the triplet form of the matrix */
	//printf ("\nA: ") ;
	//(void) umfpack_zi_report_triplet (n, n, nz, Arow, Acol, Aval, Avalz, Control) ;

	/* convert to column form */
	nz1 = max (nz,1) ;	/* ensure arrays are not of size zero. */
	Ap = new int[(n+1) * sizeof (int)];
	Ai = new int [nz1 * sizeof (int)];
	Ax = new double[nz1 * sizeof (double)];
	Az = new double[nz1 * sizeof (double)];
	if (!Ap || !Ai || !Ax || !Az) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Out of memory.}" << endl;
		goto __Exit;
	}

	/* print the triplet form of the matrix */
	//printf("\nA: ");
	//umfpack_zi_report_triplet(n, n, nz, Arow, Acol, Aval, Avalz, Control);

	status = umfpack_zi_triplet_to_col (n, n, nz, Arow, Acol, Aval, Avalz,
		Ap, Ai, Ax, Az, (int *) NULL) ;

	//Don't need any more 
	if (Arow)  { delete []Arow;  Arow  = NULL; }
	if (Acol)  { delete []Acol;  Acol  = NULL; }
	if (Aval)  { delete []Aval;  Aval  = NULL; }
	if (Avalz) { delete []Avalz; Avalz = NULL; }
	
	if (status < 0) {
		umfpack_zi_report_status (Control, status) ;
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Triplet to col failed.}" << endl;
		goto __Exit;
	}

	/* print the column-form of A */
	//printf ("\nA: ") ;
	//(void) umfpack_zi_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

	//////////////////////////////////////////////////////////////////////////

	nev = nmodes; // The number of eigenvalues to calculate
	ncv = vlen;

	Evals = new std::complex<double>[nev];
	Evecs = new std::complex<double>*[n];
	for (int i=0; i<n; i++) Evecs[i] = new std::complex<double>[nev];
	if (!Evals || !Evecs) {
		//!!Notice, be careful in editing, important for Invoker 
		cout << "[MODE SOLVER ERROR]{ARPACK:Out of memory.}" << endl;
		goto __Exit;
	}

	end = clock();
	//cout<<"ARPACK:eigs, init, time: "<<time<<" , end:"<<end<<" , end-time="<<end - time<<"ms sec="<<(end - time)/1000<<"s\n\n";
	time = end;

	//cout << "\n+++++++++++++++++++++arpack_znaupd---->begin+++++++++++++++++++" << endl;

	arpack_znaupd(Ap, Ai, Ax, Az, Control, n, nev, ncv, sigma, tol, Evals, Evecs);

	//cout << "\n+++++++++++++++++++++arpack_znaupd---->end++++++++++++++++++++" << endl;

	end = clock();
	//cout<<"ARPACK:eigs, coveraged time: "<<time<<" , end:"<<end<<" , end-time="<<end - time<<"ms sec="<<(end - time)/1000<<"s\n\n";

#if _DEBUG
	cout << "nzero: " << nz << endl;
	//cout << "shift: " << shift << endl;
#endif

	for (int i = 0; i < nev; i++) {
		d_real[i] = Evals[i]._Val[0];
		d_imag[i] = Evals[i]._Val[1];
#ifdef _DEBUG__
		cout << "eigenvalue: [" << d_real[i] << "," << d_imag[i] << "]\n";
		cout << "eigenvect[0-n]: "<< endl;
#endif		
		for (int j = 0; j < ncv; j++) { 
			v_real[i*ncv+j] = Evecs[j][i]._Val[0];
			v_imag[i*ncv+j] = Evecs[j][i]._Val[1];
#ifdef _DEBUG__
			//if (j < 10)
				cout << "[" << v_real[i*ncv+j] << "," << v_imag[i*ncv+j] <<"]\n";
#endif
			
		} 
	}
#ifdef _DEBUG__
	cout << endl;
#endif
__Exit:
	if (Ap) delete []Ap;
	if (Ai) delete []Ai;
	if (Ax) delete []Ax;
	if (Az) delete []Az;
	if (Arow) delete []Arow;
	if (Acol) delete []Acol;
	if (Aval) delete []Aval;
	if (Avalz) delete []Avalz;

	if (Evals) delete []Evals;
	if (Evecs) {
		for (int i=0; i<n; i++) { if (Evecs[i]) delete []Evecs[i]; }
		delete []Evecs;
	}	
#ifdef _DEBUG
	cout << "+++++++++++++++++++++ARPACK:eigs---->end++++++++++++++++++++" << endl;
#endif
	return 0;
}