/* ========================================================================== */
/* === umfpack_report_control =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

void umfpack_di_report_control
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_dl_report_control
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_zi_report_control
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

void umfpack_zl_report_control
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL]
) ;

/*
UMFPACK_MYFLOAT int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_di_report_control (Control) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_dl_report_control (Control) ;

complex int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_zi_report_control (Control) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;
    umfpack_zl_report_control (Control) ;

Purpose:

    Prints the current control settings.  Note that with the default print
    level, nothing is printed.  Does nothing if Control is (UMFPACK_MYFLOAT *) NULL.

Arguments:

    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;   Input argument, not modified.

	If a (UMFPACK_MYFLOAT *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_PRL]:  printing level.

	    1 or less: no output
	    2 or more: print all of Control
	    Default: 1
*/
