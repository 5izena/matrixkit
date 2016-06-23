/* ========================================================================== */
/* === umfpack_report_info ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* Copyright (c) 2005-2012 by Timothy A. Davis, http://www.suitesparse.com.   */
/* All Rights Reserved.  See ../Doc/License for License.                      */
/* -------------------------------------------------------------------------- */

void umfpack_di_report_info
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    const UMFPACK_MYFLOAT Info [UMFPACK_INFO]
) ;

void umfpack_dl_report_info
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    const UMFPACK_MYFLOAT Info [UMFPACK_INFO]
) ;

void umfpack_zi_report_info
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    const UMFPACK_MYFLOAT Info [UMFPACK_INFO]
) ;

void umfpack_zl_report_info
(
    const UMFPACK_MYFLOAT Control [UMFPACK_CONTROL],
    const UMFPACK_MYFLOAT Info [UMFPACK_INFO]
) ;

/*
UMFPACK_MYFLOAT int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    umfpack_di_report_info (Control, Info) ;

UMFPACK_MYFLOAT SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    umfpack_dl_report_info (Control, Info) ;

complex int Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    umfpack_zi_report_info (Control, Info) ;

complex SuiteSparse_long Syntax:

    #include "umfpack.h"
    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    umfpack_zl_report_info (Control, Info) ;

Purpose:

    Reports statistics from the umfpack_*_*symbolic, umfpack_*_numeric, and
    umfpack_*_*solve routines.

Arguments:

    UMFPACK_MYFLOAT Control [UMFPACK_CONTROL] ;   Input argument, not modified.

	If a (UMFPACK_MYFLOAT *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_PRL]:  printing level.

	    0 or less: no output, even when an error occurs
	    1: error messages only
	    2 or more: error messages, and print all of Info
	    Default: 1

    UMFPACK_MYFLOAT Info [UMFPACK_INFO] ;		Input argument, not modified.

	Info is an output argument of several UMFPACK routines.
	The contents of Info are printed on standard output.
*/
