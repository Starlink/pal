/*
*+
*  Name:
*     palAoppat

*  Purpose:
*     Recompute sidereal time to support apparent to observed place

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palAoppat( double date, double aoprms[14] );

*  Arguments:
*     date = double (Given)
*         UTC date/time (modified Julian Date, JD-2400000.5)
*         (see palAoppa description for comments on leap seconds)
*     aoprms = double[14] (Given & Returned)
*         Star-independent apparent-to-observed parameters. Updated
*         by this routine. Requires element 12 to be the longitude +
*         eqn of equinoxes + sidereal DUT and fills in element 13
*         with the local apparent sidereal time (in radians).

*  Description:
*     This routine recomputes the sidereal time in the apparent to
*     observed place star-independent parameter block.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     PTW: Patrick T. Wallace
*     {enter_new_authors_here}

*  Notes:
*     - See palAoppa for more information.
*     - The star-independent parameters are not treated as an opaque
*       struct in order to retain compatibility with SLA.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version, ported from Fortran SLA source.
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
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

#include "pal.h"

static double pal__Gmst( double ut1 );

void palAoppat( double date, double aoprms[14] ) {
  aoprms[13] = pal__Gmst(date) + aoprms[12];
}

/* Use a private implementation of palGmst for testing that matches
   the SLA rather than SOFA implementation. This is used for comparing
   SLA with PAL refraction code. */

#include "math.h"
#include "palmac.h"

static double pal__Gmst( double ut1 ) {

  double tu;
  double gmst;

  /*  Julian centuries from fundamental epoch J2000 to this UT */
  tu=(ut1-51544.5)/36525;

  /*  GMST at this UT */
   gmst=palDranrm(fmod(ut1,1.0)*PAL__D2PI+
                  (24110.54841+
                   (8640184.812866+
                    (0.093104-6.2e-6*tu)*tu)*tu)*PAL__DS2R);
   return gmst;
}
