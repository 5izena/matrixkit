/* ========================================================================== */
/* === umfpack_wsolve ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

int umfpack_di_wsolve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const UMFPACK_MYFLOAT Ax [ ],
    UMFPACK_MYFLOAT X [ ],
    const UMFPACK_MYFLOAT B [ ],
    void *Numeric,
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    UMFPACK_MYFLOAT Info [UMFPACK_INFO],
    int Wi [ ],
    UMFPACK_MYFLOAT W [ ]
) ;

SuiteSparse_long umfpack_dl_wsolve
(
    SuiteSparse_long sys,
    const SuiteSparse_long Ap [ ],
    const SuiteSparse_long Ai [ ],
    const UMFPACK_MYFLOAT Ax [ ],
    UMFPACK_MYFLOAT X [ ],
    const UMFPACK_MYFLOAT B [ ],
    void *Numeric,
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    UMFPACK_MYFLOAT Info [UMFPACK_INFO],
    SuiteSparse_long Wi [ ],
    UMFPACK_MYFLOAT W [ ]
) ;

int umfpack_zi_wsolve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const UMFPACK_MYFLOAT Ax [ ], const UMFPACK_MYFLOAT Az [ ],
    UMFPACK_MYFLOAT Xx [ ],	 UMFPACK_MYFLOAT Xz [ ],
    const UMFPACK_MYFLOAT Bx [ ], const UMFPACK_MYFLOAT Bz [ ],
    void *Numeric,
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    UMFPACK_MYFLOAT Info [UMFPACK_INFO],
    int Wi [ ],
    UMFPACK_MYFLOAT W [ ]
) ;

SuiteSparse_long umfpack_zl_wsolve
(
    SuiteSparse_long sys,
    const SuiteSparse_long Ap [ ],
    const SuiteSparse_long Ai [ ],
    const UMFPACK_MYFLOAT Ax [ ], const UMFPACK_MYFLOAT Az [ ],
    UMFPACK_MYFLOAT Xx [ ],	 UMFPACK_MYFLOAT Xz [ ],
    const UMFPACK_MYFLOAT Bx [ ], const UMFPACK_MYFLOAT Bz [ ],
    void *Numeric,
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    UMFPACK_MYFLOAT Info [UMFPACK_INFO],
    SuiteSparse_long Wi [ ],
    UMFPACK_MYFLOAT W [ ]
) ;

/*
UMFPACK_MYFLOAT int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, *Wi, sys ;
    UMFPACK_MYFLOAT *B, *X, *Ax, *W, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_di_wsolve (sys, Ap, Ai, Ax, X, B, Numeric,
	Control, Info, Wi, W) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    SuiteSparse_long status, *Ap, *Ai, *Wi, sys ;
    UMFPACK_MYFLOAT *B, *X, *Ax, *W, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_wsolve (sys, Ap, Ai, Ax, X, B, Numeric,
	Control, Info, Wi, W) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, *Wi, sys ;
    UMFPACK_MYFLOAT *Bx, *Bz, *Xx, *Xz, *Ax, *Az, *W,
	Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_wsolve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info, Wi, W) ;

complex SuiteSparse_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    SuiteSparse_long status, *Ap, *Ai, *Wi, sys ;
    UMFPACK_MYFLOAT *Bx, *Bz, *Xx, *Xz, *Ax, *Az, *W,
	Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_wsolve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info, Wi, W) ;

packed complex Syntax:

    Same as above, except Az, Xz, and Bz are NULL.

Purpose:

    Given LU factors computed by umfpack_*_numeric (PAQ=LU) and the
    right-hand-side, B, solve a linear system for the solution X.  Iterative
    refinement is optionally performed.  This routine is identical to
    umfpack_*_solve, except that it does not dynamically allocate any workspace.
    When you have many linear systems to solve, this routine is faster than
    umfpack_*_solve, since the workspace (Wi, W) needs to be allocated only
    once, prior to calling umfpack_*_wsolve.

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int sys ;		Input argument, not modified.
    Int Ap [n+1] ;	Input argument, not modified.
    Int Ai [nz] ;	Input argument, not modified.
    UMFPACK_MYFLOAT Ax [nz] ;	Input argument, not modified.
			Size 2*nz in packed complex case.
    UMFPACK_MYFLOAT X [n] ;	Output argument.
    UMFPACK_MYFLOAT B [n] ;	Input argument, not modified.
    void *Numeric ;	Input argument, not modified.
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;	Input argument, not modified.
    UMFPACK_MYFLOAT Info [UMFPACK_INFO] ;	Output argument.

    for complex versions:
    UMFPACK_MYFLOAT Az [nz] ;	Input argument, not modified, imaginary part
    UMFPACK_MYFLOAT Xx [n] ;	Output argument, real part.
			Size 2*n in packed complex case.
    UMFPACK_MYFLOAT Xz [n] ;	Output argument, imaginary part
    UMFPACK_MYFLOAT Bx [n] ;	Input argument, not modified, real part.
			Size 2*n in packed complex case.
    UMFPACK_MYFLOAT Bz [n] ;	Input argument, not modified, imaginary part

	The above arguments are identical to umfpack_*_solve, except that the
	error code UMFPACK_ERROR_out_of_memory will not be returned in
	Info [UMFPACK_STATUS], since umfpack_*_wsolve does not allocate any
	memory.

    Int Wi [n] ;		Workspace.
    UMFPACK_MYFLOAT W [c*n] ;		Workspace, where c is defined below.

	The Wi and W arguments are workspace used by umfpack_*_wsolve.  They
	need not be initialized on input, and their contents are undefined on
	output.  The size of W depends on whether or not iterative refinement is
	used, and which version (real or complex) is called.  Iterative
	refinement is performed if Ax=b, A'x=b, or A.'x=b is being solved,
	Control [UMFPACK_IRSTEP] > 0, and A is nonsingular.  The size of W is
	given below:

				no iter.	with iter.
				refinement	refinement
	umfpack_di_wsolve	n		5*n
	umfpack_dl_wsolve	n		5*n
	umfpack_zi_wsolve	4*n		10*n
	umfpack_zl_wsolve	4*n		10*n
*/
