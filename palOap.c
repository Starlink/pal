/*
*+
*  Name:
*     palOap

*  Purpose:
*     Observed to apparent place

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palOap ( const char *type, double ob1, double ob2, double date,
*                   double dut, double elongm, double phim, double hm,
*                   double xp, double yp, double tdk, double pmb,
*                   double rh, double wl, double tlr,
*                   double *rap, double *dap );

*  Arguments:
*     type = const char * (Given)
*        Type of coordinates - 'R', 'H' or 'A' (see below)
*     ob1 = double (Given)
*        Observed Az, HA or RA (radians; Az is N=0;E=90)
*     ob2 = double (Given)
*        Observed ZD or Dec (radians)
*     date = double (Given)
*        UTC date/time (Modified Julian Date, JD-2400000.5)
*     dut = double (Given)
*        delta UT: UT1-UTC (UTC seconds)
*     elongm = double (Given)
*        Mean longitude of the observer (radians, east +ve)
*     phim = double (Given)
*        Mean geodetic latitude of the observer (radians)
*     hm = double (Given)
*        Observer's height above sea level (metres)
*     xp = double (Given)
*        Polar motion x-coordinates (radians)
*     yp = double (Given)
*        Polar motion y-coordinates (radians)
*     tdk = double (Given)
*        Local ambient temperature (K; std=273.15)
*     pmb = double (Given)
*        Local atmospheric pressure (mb; std=1013.25)
*     rh = double (Given)
*        Local relative humidity (in the range 0.0-1.0)
*     wl = double (Given)
*        Effective wavelength (micron, e.g. 0.55)
*     tlr = double (Given)
*        Tropospheric laps rate (K/metre, e.g. 0.0065)
*     rap = double * (Given)
*        Geocentric apparent right ascension
*     dap = double * (Given)
*        Geocentric apparent declination

*  Description:
*     Observed to apparent place.

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
*     (north zero, east 90 deg) and zenith distance.  (Zenith
*     distance is used rather than elevation in order to reflect the
*     fact that no allowance is made for depression of the horizon.)
*
*     - The accuracy of the result is limited by the corrections for
*     refraction.  Providing the meteorological parameters are
*     known accurately and there are no gross local effects, the
*     predicted apparent RA,Dec should be within about 0.1 arcsec
*     for a zenith distance of less than 70 degrees.  Even at a
*     topocentric zenith distance of 90 degrees, the accuracy in
*     elevation should be better than 1 arcmin;  useful results
*     are available for a further 3 degrees, beyond which the
*     palRefro routine returns a fixed value of the refraction.
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
*     - This routine takes time to execute, due mainly to the rigorous
*     integration used to evaluate the refraction.  For processing
*     multiple stars for one location and time, call palAoppa once
*     followed by one call per star to palOapqk.  Where a range of
*     times within a limited period of a few hours is involved, and the
*     highest precision is not required, call palAoppa once, followed
*     by a call to palAoppat each time the time changes, followed by
*     one call per star to palOapqk.
*
*     - The DATE argument is UTC expressed as an MJD.  This is, strictly
*     speaking, wrong, because of leap seconds.  However, as long as
*     the delta UT and the UTC are consistent there are no
*     difficulties, except during a leap second.  In this case, the
*     start of the 61st second of the final minute should begin a new
*     MJD day and the old pre-leap delta UT should continue to be used.
*     As the 61st second completes, the MJD should revert to the start
*     of the day as, simultaneously, the delta UTC changes by one
*     second to its post-leap new value.
*
*     - The delta UT (UT1-UTC) is tabulated in IERS circulars and
*     elsewhere.  It increases by exactly one second at the end of
*     each UTC leap second, introduced in order to keep delta UT
*     within +/- 0.9 seconds.
*
*     - IMPORTANT -- TAKE CARE WITH THE LONGITUDE SIGN CONVENTION.
*     The longitude required by the present routine is east-positive,
*     in accordance with geographical convention (and right-handed).
*     In particular, note that the longitudes returned by the
*     palOBS routine are west-positive, following astronomical
*     usage, and must be reversed in sign before use in the present
*     routine.
*
*     - The polar coordinates XP,YP can be obtained from IERS
*     circulars and equivalent publications.  The maximum amplitude
*     is about 0.3 arcseconds.  If XP,YP values are unavailable,
*     use XP=YP=0D0.  See page B60 of the 1988 Astronomical Almanac
*     for a definition of the two angles.
*
*     - The height above sea level of the observing station, HM,
*     can be obtained from the Astronomical Almanac (Section J
*     in the 1988 edition), or via the routine palOBS.  If P,
*     the pressure in millibars, is available, an adequate
*     estimate of HM can be obtained from the expression
*
*            HM ~ -29.3*TSL*LOG(P/1013.25).
*
*     where TSL is the approximate sea-level air temperature in K
*     (see Astrophysical Quantities, C.W.Allen, 3rd edition,
*     section 52).  Similarly, if the pressure P is not known,
*     it can be estimated from the height of the observing
*     station, HM, as follows:
*
*            P ~ 1013.25*EXP(-HM/(29.3*TSL)).
*
*     Note, however, that the refraction is nearly proportional to the
*     pressure and that an accurate P value is important for precise
*     work.
*
*     - The azimuths etc. used by the present routine are with respect
*     to the celestial pole.  Corrections from the terrestrial pole
*     can be computed using palPolmo.

*  History:
*     2012-08-27 (TIMJ):
*        Initial version, copied from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2005 Patrick T. Wallace
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

void palOap ( const char *type, double ob1, double ob2, double date,
              double dut, double elongm, double phim, double hm,
              double xp, double yp, double tdk, double pmb,
              double rh, double wl, double tlr,
              double *rap, double *dap ) {

  double aoprms[14];

  palAoppa(date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr,
           aoprms);
  palOapqk(type,ob1,ob2,aoprms,rap,dap);

}
