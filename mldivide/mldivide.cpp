
#include <math.h>
#include <vector>
#include <time.h>
#include "mldivide.h";
#include "umfpack.h"

using namespace std;

//complex<double> csr format
struct Matrix_ZCSR 
{
	int n; //Dimension of the matrix A.
	unsigned int nnz;//the non-zero elements number
	complex<myfloat> *values;//complex array that contains the non-zero elements of a sparse matrix
	int *columns;	//Element i of the integer array columns is the number of the column that contains the i-th element in the values array.
	int *rowIndex;	//Element j of the integer array rowIndex gives the index of the element in the values array that is first non-zero element in a row j.

	Matrix_ZCSR()
	{
		n = 0;
		nnz = 0;
		values = NULL;
		columns = NULL;
		rowIndex = NULL;
	}

	Matrix_ZCSR(int n_, unsigned int nnz_, complex<myfloat> *v_, int *c_, int *r_)
	{
		n = n_;
		nnz = nnz_;
		values = v_;
		columns = c_;
		rowIndex = r_;
	}
	void operator=(const Matrix_ZCSR &mz)
	{
		n = mz.n;
		nnz = mz.nnz;
		values = mz.values;
		columns = mz.columns;
		rowIndex = mz.rowIndex;
	}

	void release()
	{
		n = 0;
		nnz = 0;
		if (values != NULL) free(values);
		if (columns != NULL) free(columns);
		if (rowIndex != NULL) free(rowIndex);
	}
} ;


/**
    Converts a sparse matrix in the COO format to the CSR format.
    nnz:(input/output)
    acoo:(input)
    rowind:(input)
    colind:(input)
    n:(input/output)
    acsr:(output)
    ja:(output)
    ia:(output)
    info:(output)
*/
void zcoocsr(int *nnz, complex<myfloat> *acoo, int *rowind, int *colind, int *n, complex<myfloat> **acsr, int **ja, int **ia, int *info)
{
    int job[6];
    {
        job[0] = 1;	//Convert from coo to CSR format
        job[1] = 1;	//One-based indexing for CSR format
        job[2] = 0;	//Zero-based indexing for coordinate format
        job[3] = 2;	//Pass the whole dense matrix to MKL routines(adns is  a whole matrix A.)
        job[4] = *nnz;//Sets number of the non-zero elements of the matrix A if job(0)=1.
        job[5] = 1;	//Only array ia is filled in for the output storage.
    }

    *ia = (int *)malloc(((*n) + 1)*sizeof(int));
    if (IS_SINGLE_PRECISION) {
        mkl_ccsrcoo(job, n, (complex<float> *)(*acsr), *ja, *ia, nnz, (complex<float> *)acoo, rowind, colind, info);
    }
    else {
        mkl_zcsrcoo(job, n, (complex<double> *)(*acsr), *ja, *ia, nnz, (complex<double> *)acoo, rowind, colind, info);
    }
    *nnz = (*ia)[(*n)] - 1;
    *ja = (int *)malloc(*nnz*sizeof(int));
    *acsr = (complex<myfloat> *)malloc((*nnz)*sizeof(complex<myfloat>));
    job[5] = 2;//It is assumed that the routine already has been called with the job(6)=1, 
    //and the user allocated the required space for storing the output arrays acsr and ja.
    if (IS_SINGLE_PRECISION) {
        mkl_ccsrcoo(job, n, (complex<float> *)(*acsr), *ja, *ia, nnz, (complex<float> *)acoo, rowind, colind, info);
    }
    else {
        mkl_zcsrcoo(job, n, (complex<double> *)(*acsr), *ja, *ia, nnz, (complex<double> *)acoo, rowind, colind, info);
    }
}

//Converts a sparse matrix in the CSR format to the COO format.
void zcsrcoo(int *n, complex<myfloat> *acsr, int *ja, int *ia, int *nnz, complex<myfloat> **acoo, int **rowind, int **colind, int *info)
{
    int job[6];
    {
        job[0] = 0;	//Convert from coo to CSR format
        job[1] = 1;	//One-based indexing for CSR format
        job[2] = 0;	//Zero-based indexing for coordinate format
        job[3] = 2;	//Pass the whole dense matrix to MKL routines(adns is  a whole matrix A.)
        job[4] = *nnz;//Sets number of the non-zero elements of the matrix A if job(0)=1.
        job[5] = 3;
    }

    *rowind = (int *)malloc((*nnz)*sizeof(int));
    *colind = (int *)malloc(*nnz*sizeof(int));
    *acoo = (complex<myfloat> *)malloc(*nnz*sizeof(complex<myfloat>));
    if (IS_SINGLE_PRECISION) {
        mkl_ccsrcoo(job, n, (complex<float> *)acsr, ja, ia, nnz, (complex<float> *)(*acoo), *rowind, *colind, info);
    }
    else {
        mkl_zcsrcoo(job, n, (complex<double> *)acsr, ja, ia, nnz, (complex<double> *)(*acoo), *rowind, *colind, info);
    }
}

//Computes product of two sparse matrices stored in the CSR format (3-array variation) with one-based indexing.
//m:Number of rows of the matrix A.
//n:Number of columns of the matrix A.
//k:Number of columns of the matrix B.
Matrix_ZCSR zcsrmultcsr(int m, int n, int k, const Matrix_ZCSR &a, const Matrix_ZCSR &b)
{

	char trans = 'N';
	int request = 0;
	int sort = 0;

	int nzmax = 0;
	complex<myfloat> *ccsr = NULL;
	int *jc = NULL;
	int *ic = NULL;
	int info = -1;

	int nnz = 0;
	ic = (int *)malloc((m+1)*sizeof(int));
	request = 1;

	if (IS_SINGLE_PRECISION) {
		mkl_ccsrmultcsr (&trans, &request, &sort, &m, &n, &k, (complex<float> *)a.values, a.columns, a.rowIndex, 
			(complex<float> *)b.values, b.columns, b.rowIndex, (complex<float> *)ccsr, jc, ic, &nzmax, &info);
	} else {
		mkl_zcsrmultcsr (&trans, &request, &sort, &m, &n, &k, (complex<double> *)a.values, a.columns, a.rowIndex, 
		(complex<double> *)b.values, b.columns, b.rowIndex, (complex<double> *)ccsr, jc, ic, &nzmax, &info);
	}

	//The value of the last element ic(m + 1) or ic(n + 1) is equal to the number of non-zero elements of the matrix C plus one.
	nnz = ic[m]-1;
	ccsr = (complex<myfloat> *)malloc(nnz*sizeof(complex<myfloat>));
	jc = (int *)malloc(nnz*sizeof(int));
	request = 2;// the routine has been called previously with the parameter request=1

	if (IS_SINGLE_PRECISION) {
		mkl_ccsrmultcsr (&trans, &request, &sort, &m, &n, &k, (complex<float> *)a.values, a.columns, a.rowIndex, 
			(complex<float> *)b.values, b.columns, b.rowIndex, (complex<float> *)ccsr, jc, ic, &nzmax, &info);
	}
	else {
		mkl_zcsrmultcsr (&trans, &request, &sort, &m, &n, &k, (complex<double> *)a.values, a.columns, a.rowIndex,
			(complex<double> *)b.values, b.columns, b.rowIndex, (complex<double> *)ccsr, jc, ic, &nzmax, &info);
	}
	return Matrix_ZCSR(m, nnz, ccsr, jc, ic);
}

//Computes the sum of two matrices stored in the CSR format (3-array variation) with one-based indexing.
//m: Number of rows of the matrix A.
//n: Number of columns of the matrix A.
Matrix_ZCSR zcsradd(int m, int n, const Matrix_ZCSR &a, const Matrix_ZCSR &b)
{
	char trans = 'N';
	int request = 0;
	complex<myfloat> beta = 1.0;
	int nzmax = 0;// The length of the arrays c and jc.This parameter is used only if request=0. 
	int sort = 3;
	complex<myfloat> *ccsr = NULL;
	int *ic = (int *)malloc((m+1)*sizeof(int));
	int *jc = NULL;
	int info = -1;

	request = 1;
	if (IS_SINGLE_PRECISION) {
		mkl_ccsradd (&trans, &request, &sort, &m, &n, (complex<float> *)a.values, a.columns, a.rowIndex, (complex<float> *)(&beta),
			(complex<float> *)b.values, b.columns, b.rowIndex, (complex<float> *)ccsr, jc, ic, &nzmax, &info);
	} else {
		mkl_zcsradd (&trans, &request, &sort, &m, &n, (complex<double> *)a.values, a.columns, a.rowIndex, (complex<double> *)(&beta),
			(complex<double> *)b.values, b.columns, b.rowIndex, (complex<double> *)ccsr, jc, ic, &nzmax, &info);
	}
	//The value of the last element ic(m + 1) or ic(n + 1) is equal to the number of non-zero elements of the matrix C plus one.
	int nnz = ic[m]-1;
	ccsr = (complex<myfloat> *)malloc(nnz*sizeof(complex<myfloat>));
	jc = (int *)malloc(nnz*sizeof(int));
	request = 2;// the routine has been called previously with the parameter request=1
	if (IS_SINGLE_PRECISION) {
		mkl_ccsradd (&trans, &request, &sort, &m, &n,  (complex<float> *)a.values, a.columns, a.rowIndex,  (complex<float> *)(&beta),
			 (complex<float> *)b.values, b.columns,	b.rowIndex,  (complex<float> *)ccsr, jc, ic, &nzmax, &info);
	} else {
		mkl_zcsradd (&trans, &request, &sort, &m, &n,  (complex<double> *)a.values, a.columns, a.rowIndex, (complex<double> *)(&beta),
			(complex<double> *)b.values, b.columns, b.rowIndex, (complex<double> *)ccsr, jc, ic, &nzmax, &info);
	}
	return Matrix_ZCSR(m, nnz, ccsr, jc, ic);
}


/* -------------------------------------------------------------------------- */
/* error: print a message and exit */
/* -------------------------------------------------------------------------- */

static void error
	(
	char *message
	)
{
	printf ("\n\n====== error: %s =====\n\n", message) ;
	exit (1) ;
}

/* -------------------------------------------------------------------------- */
/* resid: compute the residual, r = Ax-b or r = A'x=b and return maxnorm (r) */
/* A' is the complex conjugate transpose, not the array transpose */
/* -------------------------------------------------------------------------- */

static double resid
	(
	int transpose,
	int n,
	UMFPACK_MYFLOAT *b,
	UMFPACK_MYFLOAT *bz,
#ifdef ZINT
	int Ap [ ],
	int Ai [ ],
#else 
	SuiteSparse_long Ap [ ],
	SuiteSparse_long Ai [ ],
#endif
	UMFPACK_MYFLOAT Ax [ ],
	UMFPACK_MYFLOAT Az [ ],
	UMFPACK_MYFLOAT *x,
	UMFPACK_MYFLOAT *xz,
	UMFPACK_MYFLOAT **r,
	UMFPACK_MYFLOAT **rz
)
{
	int i, j, p ;
	double norm ;

	for (i = 0 ; i < n ; i++)
	{
		(*r )[i] = -b [i] ;
		(*rz)[i] = -bz[i] ;
	}
	if (transpose)
	{
		for (j = 0 ; j < n ; j++)
		{
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
				i = Ai [p] ;
				/* complex: r(j) += conj (Aij) * x (i) */
				(*r )[j] += Ax [p] * x [i] ;
				(*r )[j] += Az [p] * xz[i] ;
				(*rz)[j] -= Az [p] * x [i] ;
				(*rz)[j] += Ax [p] * xz[i] ;
			}
		}
	}
	else
	{
		for (j = 0 ; j < n ; j++)
		{
			for (p = Ap [j] ; p < Ap [j+1] ; p++)
			{
				i = Ai [p] ;
				(*r )[i] += Ax [p] * x [j] ;
				(*r )[i] -= Az [p] * xz[j] ;
				(*rz)[i] += Az [p] * x [j] ;
				(*rz)[i] += Ax [p] * xz[j] ;
			}
		}
	}
	norm = 0. ;
	for (i = 0 ; i < n ; i++)
	{
		norm = SF_MAX (SF_ABSZ ((*r )[i], (*rz )[i]), norm) ;
	}
	return (norm) ;
}

//Solves a system of linear equations for a sparse matrix in the COO format.(y = inv(A)*x)
void umfpack_sv(int n, int nz, complex<myfloat> *acoo, int *rowind, int *colind, complex<myfloat> *f, complex<myfloat> **F)
{
	//=================test==================
	
#if _DEBUG
	time_t t;
	t = clock();
	int maxitr = (int)ceil((double)2*n/1);
#endif
	//=======================================
	UMFPACK_MYFLOAT *Aval = (UMFPACK_MYFLOAT *)malloc(nz*sizeof(UMFPACK_MYFLOAT));
	UMFPACK_MYFLOAT *Avalz = (UMFPACK_MYFLOAT *)malloc(nz*sizeof(UMFPACK_MYFLOAT));
	for (int i = 0; i < nz; i++) {
		Aval[i] = (UMFPACK_MYFLOAT)(acoo[i]._Val[0]);
		Avalz[i] = (UMFPACK_MYFLOAT)(acoo[i]._Val[1]);
	}

	//int *Arow = rowind;
	//int *Acol = colind;
#ifdef ZINT
	int *Arow = (int *)malloc(nz*sizeof(int));
	int *Acol = (int *)malloc(nz*sizeof(int));
#else 
	SuiteSparse_long *Arow = (SuiteSparse_long *)malloc(nz*sizeof(SuiteSparse_long));
	SuiteSparse_long *Acol = (SuiteSparse_long *)malloc(nz*sizeof(SuiteSparse_long));
#endif
	
	for (int i = 0; i < nz; i++) {
		Arow[i] = rowind[i];
		Acol[i] = colind[i];
	}

	UMFPACK_MYFLOAT *b = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	UMFPACK_MYFLOAT *bz = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	for (int i = 0; i < n; i++) {
		b[i] = (UMFPACK_MYFLOAT)(f[i]._Val[0]);
		bz[i] = (UMFPACK_MYFLOAT)(f[i]._Val[1]);
	}

	UMFPACK_MYFLOAT *x = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	UMFPACK_MYFLOAT *xz = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	/* ---------------------------------------------------------------------- */
	/* initializations */
	/* ---------------------------------------------------------------------- */
	UMFPACK_MYFLOAT *Ax, *Az;
	UMFPACK_MYFLOAT Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], rnorm;
#ifdef ZINT
	int *Ap, *Ai;
	int nz1, status;
	int *map = NULL;
#else 
	SuiteSparse_long *Ap, *Ai;
	SuiteSparse_long nz1, status;
	SuiteSparse_long *map = NULL;
#endif
	
	void *Symbolic, *Numeric ;

#ifndef NDEBUG 
	printf ("\nUMFPACK V%d.%d (%s) demo: _zl_ version\n",
		UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;
#endif

	/* get the default control parameters */
	umfpack_zl_defaults (Control) ;

	/* change the default print level for this demo */
	/* (otherwise, nothing will print) */
	Control [UMFPACK_PRL] = 6 ;

	/* print the license agreement */
	//umfpack_zl_report_status (Control, UMFPACK_OK) ;
	//Control [UMFPACK_PRL] = 5 ;

	/* print the control parameters */
	//umfpack_zl_report_control (Control) ;


	/* ---------------------------------------------------------------------- */
	/* print A and b, and convert A to column-form */
	/* ---------------------------------------------------------------------- */

	/* print the right-hand-side */
	//printf ("\nb: ") ;
	//(void) umfpack_zl_report_vector (n, b, bz, Control) ;

	/* print the triplet form of the matrix */
	//printf ("\nA: ") ;
	//(void) umfpack_zl_report_triplet (n, n, nz, Arow, Acol, Aval, Avalz,
	//	Control) ;

	/* convert to column form */
	nz1 = SF_MAX (nz,1) ;	/* ensure arrays are not of size zero. */
#ifdef ZINT
	Ap = (int *) malloc ((n+1) * sizeof (int)) ;
	Ai = (int *) malloc (nz1 * sizeof (int)) ;
#else 
	Ap = (SuiteSparse_long *) malloc ((n+1) * sizeof (SuiteSparse_long)) ;
	Ai = (SuiteSparse_long *) malloc (nz1 * sizeof (SuiteSparse_long)) ;
#endif
	
	Ax = (UMFPACK_MYFLOAT *) malloc (nz1 * sizeof (UMFPACK_MYFLOAT)) ;
	Az = (UMFPACK_MYFLOAT *) malloc (nz1 * sizeof (UMFPACK_MYFLOAT)) ;
	if (!Ap || !Ai || !Ax || !Az)
	{
		error ("out of memory") ;
	}

	status = umfpack_zl_triplet_to_col (n, n, nz, Arow, Acol, Aval, Avalz,
		Ap, Ai, Ax, Az, map) ;

	if (status < 0)
	{
		umfpack_zl_report_status (Control, status) ;
		error ("umfpack_zl_triplet_to_col failed") ;
	}

	/* print the column-form of A */
	//printf ("\nA: ") ;
	//(void) umfpack_zl_report_matrix (n, n, Ap, Ai, Ax, Az, 1, Control) ;

	/* ---------------------------------------------------------------------- */
	/* symbolic factorization */
	/* ---------------------------------------------------------------------- */

	status = umfpack_zl_symbolic (n, n, Ap, Ai, Ax, Az, &Symbolic,
		Control, Info) ;
	if (status < 0)
	{
		umfpack_zl_report_info (Control, Info) ;
		umfpack_zl_report_status (Control, status) ;
		error ("umfpack_zl_symbolic failed") ;
	}

	/* print the symbolic factorization */
#ifndef NDEBUG 
	printf ("\nSymbolic factorization of A: \n") ;
#endif
	//(void) umfpack_zl_report_symbolic (Symbolic, Control) ;

	/* ---------------------------------------------------------------------- */
	/* numeric factorization */
	/* ---------------------------------------------------------------------- */

	status = umfpack_zl_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
		Control, Info) ;
	if (status < 0)
	{
		umfpack_zl_report_info (Control, Info) ;
		umfpack_zl_report_status (Control, status) ;
		error ("umfpack_zl_numeric failed") ;
	}

	umfpack_zl_free_symbolic (&Symbolic) ;

	/* print the numeric factorization */
#ifndef NDEBUG 
	printf ("\nNumeric factorization of A: \n") ;
#endif
	//(void) umfpack_zl_report_numeric (Numeric, Control) ;

	/* ---------------------------------------------------------------------- */
	/* solve Ax=b */
	/* ---------------------------------------------------------------------- */

	status = umfpack_zl_solve (UMFPACK_A, Ap, Ai, Ax, Az, x, xz, b, bz,
		Numeric, Control, Info) ;
	//umfpack_zl_report_info (Control, Info) ;
	//umfpack_zl_report_status (Control, status) ;
	if (status < 0)
	{
		error ("umfpack_zl_solve failed") ;
	}

	umfpack_zl_free_numeric (&Numeric) ;

	//printf ("\nx (solution of Ax=b): ") ;
	//(void) umfpack_zl_report_vector (n, x, xz, Control) ;

	for (int i = 0; i < n; i++) {
		(*F)[i]._Val[0] = x[i];
		(*F)[i]._Val[1] = xz[i];
	}

	UMFPACK_MYFLOAT *r = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	UMFPACK_MYFLOAT *rz = (UMFPACK_MYFLOAT *)malloc(n*sizeof(UMFPACK_MYFLOAT));
	rnorm = resid (0, n, b, bz, Ap, Ai, Ax, Az, x, xz, &r, &rz) ;
#ifndef NDEBUG 
	printf ("maxnorm of residual: %g\n\n", rnorm) ;
#endif
	if (r != NULL) free(r);
	if (rz != NULL) free(rz);

	if (Ap != NULL) free (Ap) ;
	if (Ai != NULL) free (Ai) ;
	if (Ax != NULL) free (Ax) ;
	if (Az != NULL) free (Az) ;

	if (x != NULL) free(x);
	if (xz != NULL) free(xz);
	if (b != NULL) free(b);
	if (bz != NULL) free(bz);
	if (Aval != NULL) free(Aval);
	if (Avalz != NULL) free(Avalz);

	if (Arow != NULL) free(Arow);
	if (Acol != NULL) free(Acol);
}

int main()
{
#if 0 //TODO: init testing data  A\f
    int Nx = 10;
    int Ny = 10;
    int m = Nx * Ny;
    int nnz = 40;

    complex<myfloat> *f = (complex<myfloat> *)malloc(Nx*Ny*sizeof(complex<myfloat>));
    complex<myfloat>  *acoo = NULL;
    int *rowind = NULL;
    int *colind = NULL;
    complex<myfloat> *F_v = (complex<myfloat> *)malloc(m*sizeof(complex<myfloat>));

    umfpack_sv(m, nnz, acoo, rowind, colind, f, &F_v);

#endif

    return 0;

}