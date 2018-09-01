/*
*+
*  Name:
*     palRanorm

*  Purpose:
*     Reduce angle to be within [0, 2*pi) (single precision)

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*      float palRanorm ( float angle );

*  Arguments:
*     angle = float (Given)
*        Angle to be reduced (radians)

*  Description:
*     The result is "angle" expressed in the range 0 to 2*pi.

*  Authors:
*     TIMJ: Tim Jenness
*     PTW: Patrick T. Wallace
*     CHJ: Christopher H. Jordan
*     {enter_new_authors_here}

*  History:
*     2018-08-27 (CHJ):
*        Initial version
*        Adapted from the Fortran SLALIB library authored by Patrick Wallace.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
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
#include "pal.h"
#include "palmac.h"

float palRanorm ( float angle ) {
  float reduced = remainderf(angle, PAL__D2PI);
  if (reduced < 0.0) reduced += PAL__D2PI;
  return reduced;
}
