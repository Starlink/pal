/*
*+
*  Name:
*     pal1.h

*  Purpose:
*     Definitions of private PAL functions

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Include file

*  Invocation:
*     #include "pal1.h"

*  Description:
*     Function prototypes for private PAL functions. Will not be
*     installed.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     2012-08-24 (TIMJ):
*        Initial version
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software; you can redistribute it and/or
*     modify it under the terms of the GNU General Public License as
*     published by the Free Software Foundation; either version 3 of
*     the License, or (at your option) any later version.
*
*     This program is distributed in the hope that it will be
*     useful, but WITHOUT ANY WARRANTY; without even the implied
*     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
*     PURPOSE. See the GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with this program; if not, write to the Free Software
*     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*     MA 02110-1301, USA.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#ifndef PAL1HDEF
#define PAL1HDEF

void pal1Atms ( double rt, double tt, double dnt, double gamal,
                double r, double * dn, double * rdndr );

void pal1Atmt ( double r0, double t0, double alpha, double gamm2,
                double delm2, double c1, double c2, double c3, double c4,
                double c5, double c6, double r,
                double *t, double *dn, double *rdndr );

#endif
