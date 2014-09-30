/*
*+
*  Name:
*     palOapqk

*  Purpose:
*     Quick observed to apparent place

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palOapqk ( const char *type, double ob1, double ob2,
*                     const  double aoprms[14], double *rap, double *dap );

*  Arguments:
*     Quick observed to apparent place.

*  Description:
*     type = const char * (Given)
*        Type of coordinates - 'R', 'H' or 'A' (see below)
*     ob1 = double (Given)
*        Observed Az, HA or RA (radians; Az is N=0;E=90)
*     ob2 = double (Given)
*        Observed ZD or Dec (radians)
*     aoprms = const double [14] (Given)
*        Star-independent apparent-to-observed parameters.
*        See palAopqk for details.
*     rap = double * (Given)
*        Geocentric apparent right ascension
*     dap = double * (Given)
*        Geocentric apparent declination

*  Authors:
*     PTW: Patrick T. Wallace
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - Only the first character of the TYPE argument is significant.
*     'R' or 'r' indicates that OBS1 and OBS2 are the observed right
*     ascension and declination;  'H' or 'h' indicates that they are
*     hour angle (west +ve) and declination;  anything else ('A' or
*     'a' is recommended) indicates that OBS1 and OBS2 are azimuth
*     (north zero, east 90 deg) and zenith distance.  (Zenith distance
*     is used rather than elevation in order to reflect the fact that
*     no allowance is made for depression of the horizon.)
*
*     - The accuracy of the result is limited by the corrections for
*     refraction.  Providing the meteorological parameters are
*     known accurately and there are no gross local effects, the
*     predicted apparent RA,Dec should be within about 0.1 arcsec
*     for a zenith distance of less than 70 degrees.  Even at a
*     topocentric zenith distance of 90 degrees, the accuracy in
*     elevation should be better than 1 arcmin;  useful results
*     are available for a further 3 degrees, beyond which the
*     palREFRO routine returns a fixed value of the refraction.
*     The complementary routines palAop (or palAopqk) and palOap
*     (or palOapqk) are self-consistent to better than 1 micro-
*     arcsecond all over the celestial sphere.
*
*     - It is advisable to take great care with units, as even
*     unlikely values of the input parameters are accepted and
*     processed in accordance with the models used.
*
*     - "Observed" Az,El means the position that would be seen by a
*     perfect theodolite located at the observer.  This is
*     related to the observed HA,Dec via the standard rotation, using
*     the geodetic latitude (corrected for polar motion), while the
*     observed HA and RA are related simply through the local
*     apparent ST.  "Observed" RA,Dec or HA,Dec thus means the
*     position that would be seen by a perfect equatorial located
*     at the observer and with its polar axis aligned to the
*     Earth's axis of rotation (n.b. not to the refracted pole).
*     By removing from the observed place the effects of
*     atmospheric refraction and diurnal aberration, the
*     geocentric apparent RA,Dec is obtained.
*
*     - Frequently, mean rather than apparent RA,Dec will be required,
*     in which case further transformations will be necessary.  The
*     palAmp etc routines will convert the apparent RA,Dec produced
*     by the present routine into an "FK5" (J2000) mean place, by
*     allowing for the Sun's gravitational lens effect, annual
*     aberration, nutation and precession.  Should "FK4" (1950)
*     coordinates be needed, the routines palFk524 etc will also
*     need to be applied.
*
*     - To convert to apparent RA,Dec the coordinates read from a
*     real telescope, corrections would have to be applied for
*     encoder zero points, gear and encoder errors, tube flexure,
*     the position of the rotator axis and the pointing axis
*     relative to it, non-perpendicularity between the mounting
*     axes, and finally for the tilt of the azimuth or polar axis
*     of the mounting (with appropriate corrections for mount
*     flexures).  Some telescopes would, of course, exhibit other
*     properties which would need to be accounted for at the
*     appropriate point in the sequence.
*
*     - The star-independent apparent-to-observed-place parameters
*     in AOPRMS may be computed by means of the palAoppa routine.
*     If nothing has changed significantly except the time, the
*     palAoppat routine may be used to perform the requisite
*     partial recomputation of AOPRMS.
*
*     - The azimuths etc used by the present routine are with respect
*     to the celestial pole.  Corrections from the terrestrial pole
*     can be computed using palPolmo.


*  History:
*     2012-08-27 (TIMJ):
*        Initial version, direct copy of Fortran SLA
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

#include "pal.h"
#include "palmac.h"

void palOapqk ( const char *type, double ob1, double ob2, const double aoprms[14],
                double *rap, double *dap ) {

  /*  breakpoint for fast/slow refraction algorithm:
   *  zd greater than arctan(4), (see palRefco routine)
   *  or vector z less than cosine(arctan(z)) = 1/sqrt(17) */
  const double zbreak = 0.242535625;

  char c;
  double c1,c2,sphi,cphi,st,ce,xaeo,yaeo,zaeo,v[3],
    xmhdo,ymhdo,zmhdo,az,sz,zdo,tz,dref,zdt,
    xaet,yaet,zaet,xmhda,ymhda,zmhda,diurab,f,hma;

  /*  coordinate type */
  c = type[0];

  /*  coordinates */
  c1 = ob1;
  c2 = ob2;

  /*  sin, cos of latitude */
  sphi = aoprms[1];
  cphi = aoprms[2];

  /*  local apparent sidereal time */
  st = aoprms[13];

  /*  standardise coordinate type */
  if (c == 'r' || c == 'R') {
    c = 'r';
  } else if (c == 'h' || c == 'H') {
    c = 'h';
  } else {
    c = 'a';
  }

  /*  if az,zd convert to cartesian (s=0,e=90) */
  if (c == 'a') {
    ce = sin(c2);
    xaeo = -cos(c1)*ce;
    yaeo = sin(c1)*ce;
    zaeo = cos(c2);
  } else {

    /*     if ra,dec convert to ha,dec */
    if (c == 'r') {
      c1 = st-c1;
    }

    /*     to cartesian -ha,dec */
    palDcs2c( -c1, c2, v );
    xmhdo = v[0];
    ymhdo = v[1];
    zmhdo = v[2];

    /*     to cartesian az,el (s=0,e=90) */
    xaeo = sphi*xmhdo-cphi*zmhdo;
    yaeo = ymhdo;
    zaeo = cphi*xmhdo+sphi*zmhdo;
  }

  /*  azimuth (s=0,e=90) */
  if (xaeo != 0.0 || yaeo != 0.0) {
    az = atan2(yaeo,xaeo);
  } else {
    az = 0.0;
  }

  /*  sine of observed zd, and observed zd */
  sz = sqrt(xaeo*xaeo+yaeo*yaeo);
  zdo = atan2(sz,zaeo);

  /*
   *  refraction
   *  ---------- */

  /*  large zenith distance? */
  if (zaeo >= zbreak) {

    /*     fast algorithm using two constant model */
    tz = sz/zaeo;
    dref = (aoprms[10]+aoprms[11]*tz*tz)*tz;

  } else {

    /*     rigorous algorithm for large zd */
    palRefro(zdo,aoprms[4],aoprms[5],aoprms[6],aoprms[7],
             aoprms[8],aoprms[0],aoprms[9],1e-8,&dref);
  }

  zdt = zdo+dref;

  /*  to cartesian az,zd */
  ce = sin(zdt);
  xaet = cos(az)*ce;
  yaet = sin(az)*ce;
  zaet = cos(zdt);

  /*  cartesian az,zd to cartesian -ha,dec */
  xmhda = sphi*xaet+cphi*zaet;
  ymhda = yaet;
  zmhda = -cphi*xaet+sphi*zaet;

  /*  diurnal aberration */
  diurab = -aoprms[3];
  f = (1.0-diurab*ymhda);
  v[0] = f*xmhda;
  v[1] = f*(ymhda+diurab);
  v[2] = f*zmhda;

  /*  to spherical -ha,dec */
  palDcc2s(v,&hma,dap);

  /*  Right Ascension */
  *rap = palDranrm(st+hma);

}
