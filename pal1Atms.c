/*
*+
*  Name:
*     pal1Atms

*  Purpose:
*     Calculate stratosphere parameters

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void pal1Atms ( double rt, double tt, double dnt, double gamal,
*                     double r, double * dn, double * rdndr );

*  Arguments:
*     rt = double (Given)
*         Height of the tropopause from centre of the Earth (metre)
*     tt = double (Given)
*         Temperature at the tropopause (K)
*     dnt = double (Given)
*         Refractive index at the tropopause
*     gamal = double (Given)
*         Constant of the atmospheric model = G*MD/R
*     r = double (Given)
*         Current distance from the centre of the Earth (metre)
*     dn = double * (Returned)
*         Refractive index at r
*     rdndr = double * (Returned)
*         r * rate the refractive index is changing at r

*  Description:
*     Refractive index and derivative with respect to height for the
*     stratosphere.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     PTW: Patrick T. Wallace
*     {enter_new_authors_here}

*  Notes:
*     - Internal routine used by palRefro.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2004 Patrick T. Wallace
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

#include <math.h>

#include "pal1.h"

void pal1Atms ( double rt, double tt, double dnt, double gamal,
                double r, double * dn, double * rdndr ) {

  double b;
  double w;

  b = gamal / tt;
  w = (dnt - 1.0) * exp( -b * (r-rt) );
  *dn = 1.0 + w;
  *rdndr = -r * b * w;

}

