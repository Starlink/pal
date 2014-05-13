/*
*+
*  Name:
*     palAopqk

*  Purpose:
*     Quick apparent to observed place

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palAopqk ( double rap, double dap, const double aoprms[14],
*                     double *aob, double *zob, double *hob,
*                     double *dob, double *rob );

*  Arguments:
*     rap = double (Given)
*        Geocentric apparent right ascension
*     dap = double (Given)
*        Geocentric apparent declination
*     aoprms = const double [14] (Given)
*        Star-independent apparent-to-observed parameters.
*
*         [0]      geodetic latitude (radians)
*         [1,2]    sine and cosine of geodetic latitude
*         [3]      magnitude of diurnal aberration vector
*         [4]      height (HM)
*         [5]      ambient temperature (T)
*         [6]      pressure (P)
*         [7]      relative humidity (RH)
*         [8]      wavelength (WL)
*         [9]      lapse rate (TLR)
*         [10,11]  refraction constants A and B (radians)
*         [12]     longitude + eqn of equinoxes + sidereal DUT (radians)
*         [13]     local apparent sidereal time (radians)
*     aob = double * (Returned)
*        Observed azimuth (radians: N=0,E=90)
*     zob = double * (Returned)
*        Observed zenith distance (radians)
*     hob = double * (Returned)
*        Observed Hour Angle (radians)
*     dob = double * (Returned)
*        Observed Declination (radians)
*     rob = double * (Returned)
*        Observed Right Ascension (radians)

*  Description:
*     Quick apparent to observed place.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - This routine returns zenith distance rather than elevation
*       in order to reflect the fact that no allowance is made for
*       depression of the horizon.
*
*     - The accuracy of the result is limited by the corrections for
*       refraction.  Providing the meteorological parameters are
*       known accurately and there are no gross local effects, the
*       observed RA,Dec predicted by this routine should be within
*       about 0.1 arcsec for a zenith distance of less than 70 degrees.
*       Even at a topocentric zenith distance of 90 degrees, the
*       accuracy in elevation should be better than 1 arcmin;  useful
*       results are available for a further 3 degrees, beyond which
*       the palRefro routine returns a fixed value of the refraction.
*       The complementary routines palAop (or palAopqk) and palOap
*       (or palOapqk) are self-consistent to better than 1 micro-
*       arcsecond all over the celestial sphere.
*
*     - It is advisable to take great care with units, as even
*       unlikely values of the input parameters are accepted and
*       processed in accordance with the models used.
*
*     - "Apparent" place means the geocentric apparent right ascension
*       and declination, which is obtained from a catalogue mean place
*       by allowing for space motion, parallax, precession, nutation,
*       annual aberration, and the Sun's gravitational lens effect.  For
*       star positions in the FK5 system (i.e. J2000), these effects can
*       be applied by means of the palMap etc routines.  Starting from
*       other mean place systems, additional transformations will be
*       needed;  for example, FK4 (i.e. B1950) mean places would first
*       have to be converted to FK5, which can be done with the
*       palFk425 etc routines.
*
*     - "Observed" Az,El means the position that would be seen by a
*       perfect theodolite located at the observer.  This is obtained
*       from the geocentric apparent RA,Dec by allowing for Earth
*       orientation and diurnal aberration, rotating from equator
*       to horizon coordinates, and then adjusting for refraction.
*       The HA,Dec is obtained by rotating back into equatorial
*       coordinates, using the geodetic latitude corrected for polar
*       motion, and is the position that would be seen by a perfect
*       equatorial located at the observer and with its polar axis
*       aligned to the Earth's axis of rotation (n.b. not to the
*       refracted pole).  Finally, the RA is obtained by subtracting
*       the HA from the local apparent ST.
*
*     - To predict the required setting of a real telescope, the
*       observed place produced by this routine would have to be
*       adjusted for the tilt of the azimuth or polar axis of the
*       mounting (with appropriate corrections for mount flexures),
*       for non-perpendicularity between the mounting axes, for the
*       position of the rotator axis and the pointing axis relative
*       to it, for tube flexure, for gear and encoder errors, and
*       finally for encoder zero points.  Some telescopes would, of
*       course, exhibit other properties which would need to be
*       accounted for at the appropriate point in the sequence.
*
*     - The star-independent apparent-to-observed-place parameters
*       in AOPRMS may be computed by means of the palAoppa routine.
*       If nothing has changed significantly except the time, the
*       palAoppat routine may be used to perform the requisite
*       partial recomputation of AOPRMS.
*
*     - At zenith distances beyond about 76 degrees, the need for
*       special care with the corrections for refraction causes a
*       marked increase in execution time.  Moreover, the effect
*       gets worse with increasing zenith distance.  Adroit
*       programming in the calling application may allow the
*       problem to be reduced.  Prepare an alternative AOPRMS array,
*       computed for zero air-pressure;  this will disable the
*       refraction corrections and cause rapid execution.  Using
*       this AOPRMS array, a preliminary call to the present routine
*       will, depending on the application, produce a rough position
*       which may be enough to establish whether the full, slow
*       calculation (using the real AOPRMS array) is worthwhile.
*       For example, there would be no need for the full calculation
*       if the preliminary call had already established that the
*       source was well below the elevation limits for a particular
*       telescope.
*
*     - The azimuths etc produced by the present routine are with
*       respect to the celestial pole.  Corrections to the terrestrial
*       pole can be computed using palPolmo.

*  History:
*     2012-08-25 (TIMJ):
*        Initial version, copied from Fortran SLA
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

#include <math.h>

#include "pal.h"

void palAopqk ( double rap, double dap, const double aoprms[14],
                double *aob, double *zob, double *hob,
                double *dob, double *rob ) {

  /*  Breakpoint for fast/slow refraction algorithm:
   *  ZD greater than arctan(4), (see palRefco routine)
   *  or vector Z less than cosine(arctan(Z)) = 1/sqrt(17) */
  const double zbreak = 0.242535625;
  int i;

  double  sphi,cphi,st,v[3],xhd,yhd,zhd,diurab,f,
    xhdt,yhdt,zhdt,xaet,yaet,zaet,azobs,
    zdt,refa,refb,zdobs,dzd,dref,ce,
    xaeo,yaeo,zaeo,hmobs,dcobs,raobs;

  /*  sin, cos of latitude */
  sphi = aoprms[1];
  cphi = aoprms[2];

  /*  local apparent sidereal time */
  st = aoprms[13];

  /*  apparent ra,dec to cartesian -ha,dec */
  palDcs2c( rap-st, dap, v );
  xhd = v[0];
  yhd = v[1];
  zhd = v[2];

  /*  diurnal aberration */
  diurab = aoprms[3];
  f = (1.0-diurab*yhd);
  xhdt = f*xhd;
  yhdt = f*(yhd+diurab);
  zhdt = f*zhd;

  /*  cartesian -ha,dec to cartesian az,el (s=0,e=90) */
  xaet = sphi*xhdt-cphi*zhdt;
  yaet = yhdt;
  zaet = cphi*xhdt+sphi*zhdt;

  /*  azimuth (n=0,e=90) */
  if (xaet == 0.0 && yaet == 0.0) {
    azobs = 0.0;
  } else {
    azobs = atan2(yaet,-xaet);
  }

  /*  topocentric zenith distance */
  zdt = atan2(sqrt(xaet*xaet+yaet*yaet),zaet);

  /*
   *  refraction
   *  ---------- */

  /*  fast algorithm using two constant model */
  refa = aoprms[10];
  refb = aoprms[11];
  palRefz(zdt,refa,refb,&zdobs);

  /*  large zenith distance? */
  if (cos(zdobs) < zbreak) {

    /*     yes: use rigorous algorithm */

    /*     initialize loop (maximum of 10 iterations) */
    i = 1;
    dzd = 1.0e1;
    while (abs(dzd) > 1e-10 && i <= 10) {

      /*        compute refraction using current estimate of observed zd */
      palRefro(zdobs,aoprms[4],aoprms[5],aoprms[6],
               aoprms[7],aoprms[8],aoprms[0],
               aoprms[9],1e-8,&dref);

      /*        remaining discrepancy */
      dzd = zdobs+dref-zdt;

      /*        update the estimate */
      zdobs = zdobs-dzd;

      /*        increment the iteration counter */
      i++;
    }
  }

  /*  to cartesian az/zd */
  ce = sin(zdobs);
  xaeo = -cos(azobs)*ce;
  yaeo = sin(azobs)*ce;
  zaeo = cos(zdobs);

  /*  cartesian az/zd to cartesian -ha,dec */
  v[0] = sphi*xaeo+cphi*zaeo;
  v[1] = yaeo;
  v[2] = -cphi*xaeo+sphi*zaeo;

  /*  to spherical -ha,dec */
  palDcc2s(v,&hmobs,&dcobs);

  /*  right ascension */
  raobs = palDranrm(st+hmobs);

  /*  return the results */
  *aob = azobs;
  *zob = zdobs;
  *hob = -hmobs;
  *dob = dcobs;
  *rob = raobs;

}
