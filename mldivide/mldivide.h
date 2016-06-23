#ifndef FDFD2D_H
#define FDFD2D_H

#include <stdint.h>
#include <complex>
using namespace std;

#define IS_SINGLE_PRECISION 1

#ifdef IS_SINGLE_PRECISION
#define myfloat float
#endif

#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))
#define SF_ABSZ(x,z) ((x) >= 0 ? (x) : -(x)) + ((z) >= 0 ? (z) : -(z))
#define SF_ABS(a)   ((a) >= 0  ? (a) : (-(a)))
#define SF_SIG(a)   ((a) >= 0  ?  1  :  -1 )
#define SF_PI (3.14159265358979323846264338328)



/**
    //inverse ERxx in coo format
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
void zcoocsr(int *nnz, complex<myfloat> *acoo, int *rowind, int *colind, int *n, complex<myfloat> **acsr, int **ja, int **ia, int *info);

void zcsrcoo(int *n, complex<myfloat> *acsr, int *ja, int *ia, int *nnz, complex<myfloat> **acoo, int **rowind, int **colind, int *info);

#ifdef __cplusplus
extern "C" {
#endif

    //The mkl_?csrmultcsr routine performs a matrix-matrix operation defined as
    //C := op(A)*B
    //op(A) is one of op(A) = A, or op(A) =A', or op(A) = conjg(A') .
    void mkl_ccsrmultcsr(char *trans, int *request, int *sort, int *m, int *n, int *k, complex<float> *a, int *ja, int *ia, complex<float> *b, int *jb, int *ib,
        complex<float> *c, int *jc, int *ic, int *nzmax, int *info);
    void mkl_zcsrmultcsr(char *trans, int *request, int *sort, int *m, int *n, int *k, complex<double> *a, int *ja, int *ia, complex<double> *b, int *jb, int *ib,
        complex<double> *c, int *jc, int *ic, int *nzmax, int *info);

    void mkl_ccsradd(char *trans, int *request, int *sort, int *m, int *n, complex<float> *a, int *ja, int *ia, complex<float> *beta,
        complex<float> *b, int *jb, int *ib, complex<float> *c, int *jc, int *ic, int *nzmax, int *info);
    void mkl_zcsradd(char *trans, int *request, int *sort, int *m, int *n, complex<double> *a, int *ja, int *ia, complex<double> *beta,
        complex<double> *b, int *jb, int *ib, complex<double> *c, int *jc, int *ic, int *nzmax, int *info);

    //The mkl_?csrsv routine solves a system of linear equations with matrix-vector operations for a sparse matrix in the CSR format:
    //y := alpha*inv(A)*x  or y := alpha*inv(A')*x,
    void mkl_ccsrsv(char *transa, int *m, complex<float> *alpha, char *matdescra, complex<float> *val, int *indx, int *pntrb, int *pntre, complex<float> *x, complex<float> *y);
    void mkl_zcsrsv(char *transa, int *m, complex<double> *alpha, char *matdescra, complex<double> *val, int *indx, int *pntrb, int *pntre, complex<double> *x, complex<double> *y);

    void mkl_ccsrcoo(int *job, int *n, complex<float> *acsr, int *ja, int *ia, int *nnz, complex<float> *acoo, int *rowind, int *colind, int *info);
    void mkl_zcsrcoo(int *job, int *n, complex<double> *acsr, int *ja, int *ia, int *nnz, complex<double> *acoo, int *rowind, int *colind, int *info);

    //The mkl_?csrgemv routine performs a matrix-vector operation defined as
    //y := A*x   or y := A'*x,
    //where:
    //x and y are vectors,
    //A is an m-by-m sparse square matrix in the CSR format (3-array variation), A' is the transpose of A.
    void mkl_ccsrgemv(char *transa, int *m, complex<float> *a, int *ia, int *ja, complex<float> *x, complex<float> *y);
    void mkl_zcsrgemv(char *transa, int *m, complex<double> *a, int *ia, int *ja, complex<double> *x, complex<double> *y);
    //A*y = x or A'*y = x
    void mkl_ccsrtrsv(char *uplo, char *transa, char *diag, int *m, complex<float> *a, int *ia, int *ja, complex<float> *x, complex<float> *y);
    void mkl_zcsrtrsv(char *uplo, char *transa, char *diag, int *m, complex<double> *a, int *ia, int *ja, complex<double> *x, complex<double> *y);
    //y := alpha*inv(A)*x or y := alpha*inv(A')*x
    void mkl_ccoosv(char *transa, int *m, complex<float> *alpha, char *matdescra, complex<float> *val, int *rowind, int *colind, int *nnz, complex<float> *x, complex<float> *y);
    void mkl_zcoosv(char *transa, int *m, complex<double> *alpha, char *matdescra, complex<double> *val, int *rowind, int *colind, int *nnz, complex<double> *x, complex<double> *y);

    //The mkl_zcsrsm routine solves a system of linear equations with matrix-matrix operations for a sparse matrix in the CSR format:
    //C := alpha*inv(A)*B or C := alpha*inv(A')*B
    void mkl_zcsrsm(char *transa, int *m, int *n, complex<double> *alpha, char *matdescra, complex<double> *val, int *indx, int *pntrb, int *pntre, complex<double> *b, int *ldb, complex<double> *c, int *ldc);


#ifdef __cplusplus
}
#endif

#endif  //FDFD2D_H