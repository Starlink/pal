/*
*+
*  Name:
*     palEcleq

*  Purpose:
*     Transform from ecliptic coordinates to J2000.0 equatorial coordinates

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palEcleq ( double dl, double db, double date,
*                     double *dr, double *dd );

*  Arguments:
*     dl = double (Given)
*        Ecliptic longitude (mean of date, IAU 1980 theory, radians)
*     db = double (Given)
*        Ecliptic latitude (mean of date, IAU 1980 theory, radians)
*     date = double (Given)
*        TT as Modified Julian Date (JD-2400000.5). The difference
*        between TT and TDB is of the order of a millisecond or two
*        (i.e. about 0.02 arc-seconds).
*     dr = double * (Returned)
*        J2000.0 mean RA (radians)
*     dd = double * (Returned)
*        J2000.0 mean Dec (Radians)

*  Description:
*     Transform from ecliptic coordinate to J2000.0 equatorial coordinates.

*  Authors:
*     PTW: Patrick T. Wallace
*     TIMJ: Tim Jenness (Cornell University)
*     {enter_new_authors_here}

*  History:
*     2014-11-18 (TIMJ):
*        Initial version
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
*     Copyright (C) 2014 Cornell University
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
#include "pal1sofa.h"

void palEcleq ( double dl, double db, double date, double *dr, double *dd ) {
  double v1[3], v2[3];
  double rmat[3][3];

  /* Spherical to Cartesian */
  eraS2c( dl, db, v1 );

  /* Ecliptic to equatorial */
  palEcmat( date, rmat );
  eraTrxp( rmat, v1, v2 );

  /* Mean of date to J2000 */
  palPrec( 2000.0, palEpj(date), rmat );
  eraTrxp( rmat, v2, v1 );

  /* Cartesian to spherical */
  eraC2s( v1, dr, dd );

  /* Express in conventional range */
  *dr = eraAnp( *dr );
  *dd = palDrange( *dd );
}
