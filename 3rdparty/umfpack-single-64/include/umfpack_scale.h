/* ========================================================================== */
/* === umfpack_scale ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

int umfpack_di_scale
(
    UMFPACK_MYFLOAT X [ ],
    const UMFPACK_MYFLOAT B [ ],
    void *Numeric
) ;

SuiteSparse_long umfpack_dl_scale
(
    UMFPACK_MYFLOAT X [ ],
    const UMFPACK_MYFLOAT B [ ],
    void *Numeric
) ;

int umfpack_zi_scale
(
    UMFPACK_MYFLOAT Xx [ ],	 UMFPACK_MYFLOAT Xz [ ],
    const UMFPACK_MYFLOAT Bx [ ], const UMFPACK_MYFLOAT Bz [ ],
    void *Numeric
) ;

SuiteSparse_long umfpack_zl_scale
(
    UMFPACK_MYFLOAT Xx [ ],	 UMFPACK_MYFLOAT Xz [ ],
    const UMFPACK_MYFLOAT Bx [ ], const UMFPACK_MYFLOAT Bz [ ],
    void *Numeric
) ;

/*
UMFPACK_MYFLOAT int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UMFPACK_MYFLOAT *B, *X ;
    status = umfpack_di_scale (X, B, Numeric) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UMFPACK_MYFLOAT *B, *X ;
    status = umfpack_dl_scale (X, B, Numeric) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UMFPACK_MYFLOAT *Bx, *Bz, *Xx, *Xz ;
    status = umfpack_zi_scale (Xx, Xz, Bx, Bz, Numeric) ;

complex SuiteSparse_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UMFPACK_MYFLOAT *Bx, *Bz, *Xx, *Xz ;
    status = umfpack_zl_scale (Xx, Xz, Bx, Bz, Numeric) ;

packed complex Syntax:

    Same as above, except both Xz and Bz are NULL.

Purpose:

    Given LU factors computed by umfpack_*_numeric (PAQ=LU, PRAQ=LU, or
    P(R\A)Q=LU), and a vector B, this routine computes X = B, X = R*B, or
    X = R\B, as appropriate.  X and B must be vectors equal in length to the
    number of rows of A.

Returns:

    The status code is returned.  UMFPACK_OK is returned if successful.
    UMFPACK_ERROR_invalid_Numeric_object is returned in the Numeric
    object is invalid.  UMFPACK_ERROR_argument_missing is returned if
    any of the input vectors are missing (X and B for the real version,
    and Xx and Bx for the complex version).

Arguments:

    UMFPACK_MYFLOAT X [n_row] ;	Output argument.
    or:
    UMFPACK_MYFLOAT Xx [n_row] ;	Output argument, real part.
			Size 2*n_row for packed complex case.
    UMFPACK_MYFLOAT Xz [n_row] ;	Output argument, imaginary part.

	The output vector X.  If either Xz or Bz are NULL, the vector
	X is in packed complex form, with the kth entry in Xx [2*k] and
	Xx [2*k+1], and likewise for B.

    UMFPACK_MYFLOAT B [n_row] ;	Input argument, not modified.
    or:
    UMFPACK_MYFLOAT Bx [n_row] ;	Input argument, not modified, real part.
			Size 2*n_row for packed complex case.
    UMFPACK_MYFLOAT Bz [n_row] ;	Input argument, not modified, imaginary part.

	The input vector B.  See above if either Xz or Bz are NULL.

    void *Numeric ;		Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.

*/
