/* ========================================================================== */
/* === umfpack_defaults ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

void umfpack_di_defaults
(
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_dl_defaults
(
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_zi_defaults
(
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_zl_defaults
(
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

/*
UMFPACK_MYFLOAT int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_di_defaults (Control) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_dl_defaults (Control) ;

complex int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_zi_defaults (Control) ;

complex SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_zl_defaults (Control) ;

Purpose:

    Sets the default control parameter settings.

Arguments:

    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;	Output argument.

	Control is set to the default control parameter settings.  You can
	then modify individual settings by changing specific entries in the
	Control array.  If Control is a (UMFPACK_MYFLOAT *) NULL pointer, then
	umfpack_*_defaults returns silently (no error is generated, since
	passing a NULL pointer for Control to any UMFPACK routine is valid).
*/
