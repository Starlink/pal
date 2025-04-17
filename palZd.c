/*
*+
*  Name:
*     palZd

*  Purpose:
*     HA, Dec to Zenith Distance

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     double palZd( double ha, double dec, double phi );

*  Arguments:
*     ha = double (Given)
*        Hour angle (radians)
*     dec = double (Given)
*        Declination (radians)
*     phi = double (Given)
*        Observatory latitude (radians)

*  Description:
*     Calculate the zenith distance for a given hour angle and declination
*     using the SOFA/ERFA routine eraHd2ae.

*  Notes:
*     - Uses eraHd2ae and converts the returned elevation to zenith
*       distance.  If elevation itself, or azimuth, is desired then
*       that routine could instead be called directly.

*  Copyright:
*     Copyright (C) 2025 East Asian Observatory.
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

#include "pal.h"
#include "palmac.h"
#include "pal1sofa.h"

double palZd( double ha, double dec, double phi ) {
  double az, el;

  eraHd2ae( ha, dec, phi, &az, &el );

  return PAL__DPIBY2 - el;
}
