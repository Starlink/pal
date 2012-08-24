/*
*+
*  Name:
*     pal1Atmt

*  Purpose:
*     Calculate troposphere parameters

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void pal1Atmt ( double r0, double t0, double alpha, double gamm2,
*                     double delm2, double c1, double c2, double c3,
*                     double c4, double c5, double c6, double r,
*                     double *t, double *dn, double *rdndr );

*  Arguments:
*     r0 = double (Given)
*         Height of observer from centre of the Earth (metre)
*     t0 = double (Given)
*         Temperature of the observer (K)
*     alpha = double (Given)
*         Alpha (see HMNAO paper)
*     gamm2 = double (Given)
*         Gamma minus 2 (see HMNAO paper)
*     delm2 = double (Given)
*         Delta minus 2 (see HMNAO paper)
*     c1 = double (Given)
*         Useful term (see palRefro source)
*     c2 = double (Given)
*         Useful term (see palRefro source)
*     c3 = double (Given)
*         Useful term (see palRefro source)
*     c4 = double (Given)
*         Useful term (see palRefro source)
*     c5 = double (Given)
*         Useful term (see palRefro source)
*     c6 = double (Given)
*         Useful term (see palRefro source)
*     r = double (Given)
*         Current distance from the centre of the Earth (metre)
*     t = double * (Returned)
*         Temperature at r (K)
*     dn = double * (Returned)
*         Refractive index at r.
*     rdndr = double * (Returned)
*         r * rate the refractive index is changing at r.

*  Description:
*     Refractive index and derivative with respect to height for
*     the troposphere.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - Internal routine used by palRefro
*     - Note that in the optical case c5 and c6 are zero.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version, copied from Fortran SLA source.
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

#include <math.h>

#include "palmac.h"
#include "pal1.h"

void pal1Atmt ( double r0, double t0, double alpha, double gamm2,
                double delm2, double c1, double c2, double c3, double c4,
                double c5, double c6, double r,
                double *t, double *dn, double *rdndr ) {

  double tt0;
  double tt0gm2;
  double tt0dm2;

  *t = DMAX( DMIN( t0 - alpha*(r-r0), 320.0), 100.0 );
  tt0 = *t / t0;
  tt0gm2 = pow( tt0, gamm2 );
  tt0dm2 = pow( tt0, delm2 );
  *dn = 1.0 + ( c1 * tt0gm2 - ( c2 - c5 / *t ) * tt0dm2 ) * tt0;
  *rdndr = r * ( -c3 * tt0gm2 + ( c4 - c6 / tt0 ) * tt0dm2 );

}
