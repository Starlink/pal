/*
*+
*  Name:
*     palTp2s

*  Purpose:
*     Tangent plane to spherical coordinates

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDtp2s( float xi, float eta, float raz, float decz,
*               float *ra, float *dec);

*  Arguments:
*     xi = float (Given)
*        First rectangular coordinate on tangent plane (radians)
*     eta = float (Given)
*        Second rectangular coordinate on tangent plane (radians)
*     raz = float (Given)
*        RA spherical coordinate of tangent point (radians)
*     decz = float (Given)
*        Dec spherical coordinate of tangent point (radians)
*     ra = float * (Returned)
*        RA spherical coordinate of point to be projected (radians)
*     dec = float * (Returned)
*        Dec spherical coordinate of point to be projected (radians)

*  Description:
*     Transform tangent plane coordinates into spherical.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     CHJ: Christopher H. Jordan
*     {enter_new_authors_here}

*  History:
*     2019-05-29 (CHJ):
*        Initial version, derived from palDtp2s.c
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#include "pal.h"
#include "pal1sofa.h"

#include <math.h>

void
palTp2s ( float xi, float eta, float raz, float decz,
           float *ra, float *dec ) {

  float cdecz;
  float denom;
  float sdecz;
  float d;

  sdecz = sinf(decz);
  cdecz = cosf(decz);
  denom = cdecz - eta * sdecz;
  d = atan2f(xi, denom) + raz;
  *ra = palRanorm(d);
  *dec = atan2f(sdecz + eta * cdecz, sqrt(xi * xi + denom * denom));

  return;
}
