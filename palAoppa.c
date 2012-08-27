/*
*+
*  Name:
*     palAoppa

*  Purpose:
*     Precompute apparent to observed place parameters

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palAoppa ( double date, double dut, double elongm, double phim,
*                     double hm, double xp, double yp, double tdk, double pmb,
*                     double rh, double wl, double tlr, double aoprms[14] );

*  Arguments:
*     date = double (Given)
*        UTC date/time (modified Julian Date, JD-2400000.5)
*     dut = double (Given)
*        delta UT:  UT1-UTC (UTC seconds)
*     elongm = double (Given)
*        mean longitude of the observer (radians, east +ve)
*     phim = double (Given)
*        mean geodetic latitude of the observer (radians)
*     hm = double (Given)
*        observer's height above sea level (metres)
*     xp = double (Given)
*        polar motion x-coordinate (radians)
*     yp = double (Given)
*        polar motion y-coordinate (radians)
*     tdk = double (Given)
*        local ambient temperature (K; std=273.15)
*     pmb = double (Given)
*        local atmospheric pressure (mb; std=1013.25)
*     rh = double (Given)
*        local relative humidity (in the range 0.0-1.0)
*     wl = double (Given)
*        effective wavelength (micron, e.g. 0.55)
*     tlr = double (Given)
*        tropospheric lapse rate (K/metre, e.g. 0.0065)
*     aoprms = double [14] (Returned)
*        Star-independent apparent-to-observed parameters
*
*         (0)      geodetic latitude (radians)
*         (1,2)    sine and cosine of geodetic latitude
*         (3)      magnitude of diurnal aberration vector
*         (4)      height (hm)
*         (5)      ambient temperature (tdk)
*         (6)      pressure (pmb)
*         (7)      relative humidity (rh)
*         (8)      wavelength (wl)
*         (9)     lapse rate (tlr)
*         (10,11)  refraction constants A and B (radians)
*         (12)     longitude + eqn of equinoxes + sidereal DUT (radians)
*         (13)     local apparent sidereal time (radians)

*  Description:
*     Precompute apparent to observed place parameters required by palAopqk
*     and palOapqk.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - It is advisable to take great care with units, as even
*       unlikely values of the input parameters are accepted and
*       processed in accordance with the models used.
*
*     - The DATE argument is UTC expressed as an MJD.  This is,
*       strictly speaking, improper, because of leap seconds.  However,
*       as long as the delta UT and the UTC are consistent there
*       are no difficulties, except during a leap second.  In this
*       case, the start of the 61st second of the final minute should
*       begin a new MJD day and the old pre-leap delta UT should
*       continue to be used.  As the 61st second completes, the MJD
*       should revert to the start of the day as, simultaneously,
*       the delta UTC changes by one second to its post-leap new value.
*
*     - The delta UT (UT1-UTC) is tabulated in IERS circulars and
*       elsewhere.  It increases by exactly one second at the end of
*       each UTC leap second, introduced in order to keep delta UT
*       within +/- 0.9 seconds.
*
*     - IMPORTANT -- TAKE CARE WITH THE LONGITUDE SIGN CONVENTION.
*       The longitude required by the present routine is east-positive,
*       in accordance with geographical convention (and right-handed).
*       In particular, note that the longitudes returned by the
*       palObs routine are west-positive, following astronomical
*       usage, and must be reversed in sign before use in the present
*       routine.
*
*     - The polar coordinates XP,YP can be obtained from IERS
*       circulars and equivalent publications.  The maximum amplitude
*       is about 0.3 arcseconds.  If XP,YP values are unavailable,
*       use XP=YP=0.0.  See page B60 of the 1988 Astronomical Almanac
*       for a definition of the two angles.
*
*     - The height above sea level of the observing station, HM,
*       can be obtained from the Astronomical Almanac (Section J
*       in the 1988 edition), or via the routine palObs.  If P,
*       the pressure in millibars, is available, an adequate
*       estimate of HM can be obtained from the expression
*
*             HM ~ -29.3*TSL*log(P/1013.25).
*
*       where TSL is the approximate sea-level air temperature in K
*       (see Astrophysical Quantities, C.W.Allen, 3rd edition,
*       section 52).  Similarly, if the pressure P is not known,
*       it can be estimated from the height of the observing
*       station, HM, as follows:
*
*             P ~ 1013.25*exp(-HM/(29.3*TSL)).
*
*       Note, however, that the refraction is nearly proportional to the
*       pressure and that an accurate P value is important for precise
*       work.
*
*     - Repeated, computationally-expensive, calls to palAoppa for
*       times that are very close together can be avoided by calling
*       palAoppa just once and then using palAoppat for the subsequent
*       times.  Fresh calls to palAoppa will be needed only when
*       changes in the precession have grown to unacceptable levels or
*       when anything affecting the refraction has changed.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version, ported directly from Fortran SLA.
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

#include "math.h"

#include "pal.h"
#include "palmac.h"

/* These are local SLA implementations to aid in testing. Switch
 * to native PAL implementations when tests are complete. */
static void pal__Geoc( double p, double h, double *r, double * z );
static void pal__Nutc ( double date, double * dpsi, double *deps, double * eps0 );
static double pal__Eqeqx( double date );

void palAoppa ( double date, double dut, double elongm, double phim,
                double hm, double xp, double yp, double tdk, double pmb,
                double rh, double wl, double tlr, double aoprms[14] ) {

  /* Constants */
  const double C = 173.14463331; /* Speed of light (AU per day) */
  const double SOLSID = 1.0027379093; /* Ratio between solar and sidereal time */

  /* Local variables */
  double cphim,xt,yt,zt,xc,yc,zc,elong,phi,uau,vau;

  /*  Observer's location corrected for polar motion */
  cphim = cos(phim);
  xt = cos(elongm)*cphim;
  yt = sin(elongm)*cphim;
  zt = sin(phim);
  xc = xt-xp*zt;
  yc = yt+yp*zt;
  zc = xp*xt-yp*yt+zt;
  if (xc == 0.0 && yc == 0.0) {
    elong = 0.0;
  } else {
    elong = atan2(yc,xc);
  }
  phi = atan2(zc,sqrt(xc*xc+yc*yc));
  aoprms[0] = phi;
  aoprms[1] = sin(phi);
  aoprms[2] = cos(phi);

  /*  magnitude of the diurnal aberration vector */
  pal__Geoc(phi,hm,&uau,&vau);
  aoprms[3] = PAL__D2PI*uau*SOLSID/C;

  /*  copy the refraction parameters and compute the a & b constants */
  aoprms[4] = hm;
  aoprms[5] = tdk;
  aoprms[6] = pmb;
  aoprms[7] = rh;
  aoprms[8] = wl;
  aoprms[9] = tlr;
  palRefco(hm,tdk,pmb,rh,wl,phi,tlr,1e-10,
           &aoprms[10],&aoprms[11]);

  /*  longitude + equation of the equinoxes + sidereal equivalent of DUT
   *  (ignoring change in equation of the equinoxes between UTC and TDB) */
  aoprms[12] = elong+pal__Eqeqx(date)+dut*SOLSID*PAL__DS2R;

  /*  sidereal time */
  palAoppat(date,aoprms);

}

/* Private reimplementation of slaEqeqx for testing the algorithm */

#include <math.h>

static void pal__Geoc( double p, double h, double *r, double * z ) {
  /*  earth equatorial radius (metres) */
  const double A0=6378140.0;

  /*  reference spheroid flattening factor and useful function */
  const double f = 1.0/298.257;
  double b;

  /*  astronomical unit in metres */
  const double AU = 1.49597870e11;

  double sp,cp,c,s;

  b = pow( 1.0-f, 2.0 );

  /*  geodetic to geocentric conversion */
  sp = sin(p);
  cp = cos(p);
  c = 1.0/sqrt(cp*cp+b*sp*sp);
  s = b*c;
  *r = (A0*c+h)*cp/AU;
  *z = (A0*s+h)*sp/AU;

}

static double pal__Eqeqx( double date ) {

  const double T2AS=1296000.0;

  double sla_eqeqx;
  double t, om, dpsi, deps, eps0;

  /*  interval between basic epoch j2000.0 and current epoch (jc) */
  t=(date-51544.5)/36525.0;

  /*  longitude of the mean ascending node of the lunar orbit on the
   *   ecliptic, measured from the mean equinox of date */
  om=PAL__DAS2R*(450160.280+(-5.0*T2AS-482890.539
                               +(7.455+0.008*t)*t)*t);

  /*  nutation */
  pal__Nutc(date,&dpsi,&deps,&eps0);

  /*  equation of the equinoxes */
  sla_eqeqx=dpsi*cos(eps0)+PAL__DAS2R*(0.00264*sin(om)+
                                 0.000063*sin(om+om));

  return sla_eqeqx;
}

#include "palmac.h"

static void pal__Nutc ( double date, double * dpsi, double *deps, double * eps0 ) {

  const double DJC = 36525.0;
  const double DJM0 = 51544.5;
  const double TURNAS = 1296000.0;

  #define NTERMS 194

  int i, j;
  double t,el,elp,f,d,om,ve,ma,ju,sa,theta,c,s,dp,de;

  int na[         194 ][9] = {
    {            0 ,            0 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            0 ,           -2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            2 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,           -2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,           -1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            3 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,           -1 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,           -2 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            1 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            1 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            1 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            1 ,            0 ,            1 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            4 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -2 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            0 ,           -2 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            4 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            4 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,           -1 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            3 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            3 ,            0 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            4 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            4 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            3 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            0 ,            4 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,           -1 ,            0 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -1 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            1 ,            0 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            4 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,           -1 ,            0 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            0 ,           -2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            1 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            3 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,           -2 ,            4 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,           -1 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,           -1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            3 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            1 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            1 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            2 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,           -1 ,            0 ,            2 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,           -2 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            1 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            4 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,           -1 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,            1 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            2 ,           -1 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            2 ,           -2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            1 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,            4 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            2 ,           -2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            1 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            2 ,           -1 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,           -1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            4 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            1 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            2 ,            1 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            4 ,           -2 ,            2 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            1 ,            0 ,            0 ,           -1 ,            0 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            4 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            1 ,            0 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,            0 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            4 ,            0 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {           -1 ,            0 ,            0 ,            4 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            0 ,            2 ,            2 ,            1 ,            0 ,            0 ,            0 ,            0  },
    {            2 ,            1 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            5 ,           -5 ,            5 ,           -3 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            2 ,            0  },
    {            0 ,            0 ,            1 ,           -1 ,            1 ,            0 ,            0 ,           -1 ,            0  },
    {            0 ,            0 ,           -1 ,            1 ,           -1 ,            1 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -1 ,            1 ,            0 ,            0 ,            2 ,            0 ,            0  },
    {            0 ,            0 ,            3 ,           -3 ,            3 ,            0 ,            0 ,           -1 ,            0  },
    {            0 ,            0 ,           -8 ,            8 ,           -7 ,            5 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -1 ,            1 ,           -1 ,            0 ,            2 ,            0 ,            0  },
    {            0 ,            0 ,           -2 ,            2 ,           -2 ,            2 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -6 ,            6 ,           -6 ,            4 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -2 ,            2 ,           -2 ,            0 ,            8 ,           -3 ,            0  },
    {            0 ,            0 ,            6 ,           -6 ,            6 ,            0 ,           -8 ,            3 ,            0  },
    {            0 ,            0 ,            4 ,           -4 ,            4 ,           -2 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -3 ,            3 ,           -3 ,            2 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            4 ,           -4 ,            3 ,            0 ,           -8 ,            3 ,            0  },
    {            0 ,            0 ,           -4 ,            4 ,           -5 ,            0 ,            8 ,           -3 ,            0  },
    {            0 ,            0 ,            0 ,            0 ,            0 ,            2 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,           -4 ,            4 ,           -4 ,            3 ,            0 ,            0 ,            0  },
    {            0 ,            1 ,           -1 ,            1 ,           -1 ,            0 ,            0 ,            1 ,            0  },
    {            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            1 ,            0  },
    {            0 ,            0 ,            1 ,           -1 ,            1 ,            1 ,            0 ,            0 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            2 ,            0 ,           -2 ,            0 ,            0  },
    {            0 ,           -1 ,           -7 ,            7 ,           -7 ,            5 ,            0 ,            0 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,            2 ,            0 ,            0 ,           -2 ,            0  },
    {           -2 ,            0 ,            2 ,            0 ,            1 ,            0 ,            0 ,           -3 ,            0  },
    {            0 ,            0 ,            2 ,           -2 ,            2 ,            0 ,            0 ,           -2 ,            0  },
    {            0 ,            0 ,            1 ,           -1 ,            1 ,            0 ,            0 ,            1 ,            0  },
    {            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            2  },
    {            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            0 ,            1  },
    {            2 ,            0 ,           -2 ,            0 ,           -2 ,            0 ,            0 ,            3 ,            0  },
    {            0 ,            0 ,            1 ,           -1 ,            1 ,            0 ,            0 ,           -2 ,            0  },
    {            0 ,            0 ,           -7 ,            7 ,           -7 ,            5 ,            0 ,            0 ,            0  }
  };
  double psi[         194 ][4] = {
    {    3341.5000000000000      ,    17206241.800000001      ,    3.1000000000000001      ,    17409.500000000000       },
    {   -1716.8000000000000      ,   -1317185.3000000000      ,    1.3999999999999999      ,   -156.80000000000001       },
    {    285.69999999999999      ,   -227667.00000000000      ,   0.29999999999999999      ,   -23.500000000000000       },
    {   -68.599999999999994      ,   -207448.00000000000      ,    0.0000000000000000      ,   -21.399999999999999       },
    {    950.29999999999995      ,    147607.89999999999      ,   -2.2999999999999998      ,   -355.00000000000000       },
    {   -66.700000000000003      ,   -51689.099999999999      ,   0.20000000000000001      ,    122.59999999999999       },
    {   -108.59999999999999      ,    71117.600000000006      ,    0.0000000000000000      ,    7.0000000000000000       },
    {    35.600000000000001      ,   -38740.199999999997      ,   0.10000000000000001      ,   -36.200000000000003       },
    {    85.400000000000006      ,   -30127.599999999999      ,    0.0000000000000000      ,   -3.1000000000000001       },
    {    9.0000000000000000      ,    21583.000000000000      ,   0.10000000000000001      ,   -50.299999999999997       },
    {    22.100000000000001      ,    12822.799999999999      ,    0.0000000000000000      ,    13.300000000000001       },
    {    3.3999999999999999      ,    12350.799999999999      ,    0.0000000000000000      ,    1.3000000000000000       },
    {   -21.100000000000001      ,    15699.400000000000      ,    0.0000000000000000      ,    1.6000000000000001       },
    {    4.2000000000000002      ,    6313.8000000000002      ,    0.0000000000000000      ,    6.2000000000000002       },
    {   -22.800000000000001      ,    5796.8999999999996      ,    0.0000000000000000      ,    6.0999999999999996       },
    {    15.699999999999999      ,   -5961.1000000000004      ,    0.0000000000000000      ,  -0.59999999999999998       },
    {    13.100000000000000      ,   -5159.1000000000004      ,    0.0000000000000000      ,   -4.5999999999999996       },
    {    1.8000000000000000      ,    4592.6999999999998      ,    0.0000000000000000      ,    4.5000000000000000       },
    {   -17.500000000000000      ,    6336.0000000000000      ,    0.0000000000000000      ,   0.69999999999999996       },
    {    16.300000000000001      ,   -3851.0999999999999      ,    0.0000000000000000      ,  -0.40000000000000002       },
    {   -2.7999999999999998      ,    4771.6999999999998      ,    0.0000000000000000      ,   0.50000000000000000       },
    {    13.800000000000001      ,   -3099.3000000000002      ,    0.0000000000000000      ,  -0.29999999999999999       },
    {   0.20000000000000001      ,    2860.3000000000002      ,    0.0000000000000000      ,   0.29999999999999999       },
    {    1.3999999999999999      ,    2045.3000000000000      ,    0.0000000000000000      ,    2.0000000000000000       },
    {   -8.5999999999999996      ,    2922.5999999999999      ,    0.0000000000000000      ,   0.29999999999999999       },
    {   -7.7000000000000002      ,    2587.9000000000001      ,    0.0000000000000000      ,   0.20000000000000001       },
    {    8.8000000000000007      ,   -1408.0999999999999      ,    0.0000000000000000      ,    3.7000000000000002       },
    {    1.3999999999999999      ,    1517.5000000000000      ,    0.0000000000000000      ,    1.5000000000000000       },
    {   -1.8999999999999999      ,   -1579.7000000000000      ,    0.0000000000000000      ,    7.7000000000000002       },
    {    1.3000000000000000      ,   -2178.5999999999999      ,    0.0000000000000000      ,  -0.20000000000000001       },
    {   -4.7999999999999998      ,    1286.8000000000000      ,    0.0000000000000000      ,    1.3000000000000000       },
    {    6.2999999999999998      ,    1267.2000000000000      ,    0.0000000000000000      ,   -4.0000000000000000       },
    {   -1.0000000000000000      ,    1669.3000000000000      ,    0.0000000000000000      ,   -8.3000000000000007       },
    {    2.3999999999999999      ,   -1020.0000000000000      ,    0.0000000000000000      ,  -0.90000000000000002       },
    {    4.5000000000000000      ,   -766.89999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.1000000000000001      ,    756.50000000000000      ,    0.0000000000000000      ,   -1.7000000000000000       },
    {   -1.3999999999999999      ,   -1097.3000000000000      ,    0.0000000000000000      ,  -0.50000000000000000       },
    {    2.6000000000000001      ,   -663.00000000000000      ,    0.0000000000000000      ,  -0.59999999999999998       },
    {   0.80000000000000004      ,   -714.10000000000002      ,    0.0000000000000000      ,    1.6000000000000001       },
    {   0.40000000000000002      ,   -629.89999999999998      ,    0.0000000000000000      ,  -0.59999999999999998       },
    {   0.29999999999999999      ,    580.39999999999998      ,    0.0000000000000000      ,   0.59999999999999998       },
    {   -1.6000000000000001      ,    577.29999999999995      ,    0.0000000000000000      ,   0.50000000000000000       },
    {  -0.90000000000000002      ,    644.39999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    2.2000000000000002      ,   -534.00000000000000      ,    0.0000000000000000      ,  -0.50000000000000000       },
    {   -2.5000000000000000      ,    493.30000000000001      ,    0.0000000000000000      ,   0.50000000000000000       },
    {  -0.10000000000000001      ,   -477.30000000000001      ,    0.0000000000000000      ,   -2.3999999999999999       },
    {  -0.90000000000000002      ,    735.00000000000000      ,    0.0000000000000000      ,   -1.7000000000000000       },
    {   0.69999999999999996      ,    406.19999999999999      ,    0.0000000000000000      ,   0.40000000000000002       },
    {   -2.7999999999999998      ,    656.89999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.59999999999999998      ,    358.00000000000000      ,    0.0000000000000000      ,    2.0000000000000000       },
    {  -0.69999999999999996      ,    472.50000000000000      ,    0.0000000000000000      ,   -1.1000000000000001       },
    {  -0.10000000000000001      ,   -300.50000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.2000000000000000      ,    435.10000000000002      ,    0.0000000000000000      ,   -1.0000000000000000       },
    {    1.8000000000000000      ,   -289.39999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.59999999999999998      ,   -422.60000000000002      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.80000000000000004      ,   -287.60000000000002      ,    0.0000000000000000      ,   0.59999999999999998       },
    {   -38.600000000000001      ,   -392.30000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.69999999999999996      ,   -281.80000000000001      ,    0.0000000000000000      ,   0.59999999999999998       },
    {   0.59999999999999998      ,   -405.69999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.2000000000000000      ,    229.00000000000000      ,    0.0000000000000000      ,   0.20000000000000001       },
    {    1.1000000000000001      ,   -264.30000000000001      ,    0.0000000000000000      ,   0.50000000000000000       },
    {  -0.69999999999999996      ,    247.90000000000001      ,    0.0000000000000000      ,  -0.50000000000000000       },
    {  -0.20000000000000001      ,    218.00000000000000      ,    0.0000000000000000      ,   0.20000000000000001       },
    {   0.59999999999999998      ,   -339.00000000000000      ,    0.0000000000000000      ,   0.80000000000000004       },
    {  -0.69999999999999996      ,    198.69999999999999      ,    0.0000000000000000      ,   0.20000000000000001       },
    {   -1.5000000000000000      ,    334.00000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,    334.00000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,   -198.09999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -106.59999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.50000000000000000      ,    165.80000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    134.80000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.90000000000000002      ,   -151.59999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -129.69999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.80000000000000004      ,   -132.80000000000001      ,    0.0000000000000000      ,  -0.10000000000000001       },
    {   0.50000000000000000      ,   -140.69999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    138.40000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    129.00000000000000      ,    0.0000000000000000      ,  -0.29999999999999999       },
    {   0.50000000000000000      ,   -121.20000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.29999999999999999      ,    114.50000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    101.80000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -3.6000000000000001      ,   -101.90000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.80000000000000004      ,   -109.40000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.20000000000000001      ,   -97.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.69999999999999996      ,    157.30000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.20000000000000001      ,   -83.299999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.29999999999999999      ,    93.299999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    92.099999999999994      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.50000000000000000      ,    133.59999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    81.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    123.90000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.29999999999999999      ,    128.09999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,    74.099999999999994      ,    0.0000000000000000      ,  -0.29999999999999999       },
    {  -0.20000000000000001      ,   -70.299999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.40000000000000002      ,    66.599999999999994      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -66.700000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.69999999999999996      ,    69.299999999999997      ,    0.0000000000000000      ,  -0.29999999999999999       },
    {    0.0000000000000000      ,   -70.400000000000006      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    101.50000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.50000000000000000      ,   -69.099999999999994      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    58.500000000000000      ,    0.0000000000000000      ,   0.20000000000000001       },
    {   0.10000000000000001      ,   -94.900000000000006      ,    0.0000000000000000      ,   0.20000000000000001       },
    {    0.0000000000000000      ,    52.899999999999999      ,    0.0000000000000000      ,  -0.20000000000000001       },
    {   0.10000000000000001      ,    86.700000000000003      ,    0.0000000000000000      ,  -0.20000000000000001       },
    {  -0.10000000000000001      ,   -59.200000000000003      ,    0.0000000000000000      ,   0.20000000000000001       },
    {   0.29999999999999999      ,   -58.799999999999997      ,    0.0000000000000000      ,   0.10000000000000001       },
    {  -0.29999999999999999      ,    49.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    56.899999999999999      ,    0.0000000000000000      ,  -0.10000000000000001       },
    {   0.29999999999999999      ,   -50.200000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    53.399999999999999      ,    0.0000000000000000      ,  -0.10000000000000001       },
    {   0.10000000000000001      ,   -76.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    45.299999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -46.799999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.20000000000000001      ,   -44.600000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.20000000000000001      ,   -48.700000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -46.799999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -42.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    46.399999999999999      ,    0.0000000000000000      ,  -0.10000000000000001       },
    {   0.20000000000000001      ,   -67.299999999999997      ,    0.0000000000000000      ,   0.10000000000000001       },
    {    0.0000000000000000      ,   -65.799999999999997      ,    0.0000000000000000      ,   0.20000000000000001       },
    {  -0.10000000000000001      ,   -43.899999999999999      ,    0.0000000000000000      ,   0.29999999999999999       },
    {    0.0000000000000000      ,   -38.899999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.29999999999999999      ,    63.899999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    41.200000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -36.100000000000001      ,    0.0000000000000000      ,   0.20000000000000001       },
    {  -0.29999999999999999      ,    58.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    36.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -39.700000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -57.700000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    33.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    36.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    55.700000000000003      ,    0.0000000000000000      ,  -0.10000000000000001       },
    {   0.10000000000000001      ,   -35.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -31.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    30.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.29999999999999999      ,    49.200000000000003      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    49.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    33.600000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -33.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -31.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    28.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -25.199999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -26.199999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    41.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    24.500000000000000      ,    0.0000000000000000      ,   0.10000000000000001       },
    {   -16.199999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -22.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    23.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    37.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.20000000000000001      ,   -25.699999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    25.199999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -24.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    24.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -20.699999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -20.800000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,    33.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    32.899999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -32.600000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    19.899999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.10000000000000001      ,    19.600000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -18.699999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -19.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.10000000000000001      ,   -28.600000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    4.0000000000000000      ,    178.80000000000001      ,   -11.800000000000001      ,   0.29999999999999999       },
    {    39.799999999999997      ,   -107.30000000000000      ,   -5.5999999999999996      ,   -1.0000000000000000       },
    {    9.9000000000000004      ,    164.00000000000000      ,   -4.0999999999999996      ,   0.10000000000000001       },
    {   -4.7999999999999998      ,   -135.30000000000001      ,   -3.3999999999999999      ,  -0.10000000000000001       },
    {    50.500000000000000      ,    75.000000000000000      ,    1.3999999999999999      ,   -1.2000000000000000       },
    {   -1.1000000000000001      ,   -53.500000000000000      ,    1.3000000000000000      ,    0.0000000000000000       },
    {   -45.000000000000000      ,   -2.3999999999999999      ,  -0.40000000000000002      ,    6.5999999999999996       },
    {   -11.500000000000000      ,   -61.000000000000000      ,  -0.90000000000000002      ,   0.40000000000000002       },
    {    4.4000000000000004      ,   -68.400000000000006      ,   -3.3999999999999999      ,    0.0000000000000000       },
    {    7.7000000000000002      ,   -47.100000000000001      ,   -4.7000000000000002      ,   -1.0000000000000000       },
    {   -42.899999999999999      ,   -12.600000000000000      ,   -1.2000000000000000      ,    4.2000000000000002       },
    {   -42.799999999999997      ,    12.699999999999999      ,   -1.2000000000000000      ,   -4.2000000000000002       },
    {   -7.5999999999999996      ,   -44.100000000000001      ,    2.1000000000000001      ,  -0.50000000000000000       },
    {   -64.099999999999994      ,    1.7000000000000000      ,   0.20000000000000001      ,    4.5000000000000000       },
    {    36.399999999999999      ,   -10.400000000000000      ,    1.0000000000000000      ,    3.5000000000000000       },
    {    35.600000000000001      ,    10.199999999999999      ,    1.0000000000000000      ,   -3.5000000000000000       },
    {   -1.7000000000000000      ,    39.500000000000000      ,    2.0000000000000000      ,    0.0000000000000000       },
    {    50.899999999999999      ,   -8.1999999999999993      ,  -0.80000000000000004      ,   -5.0000000000000000       },
    {    0.0000000000000000      ,    52.299999999999997      ,    1.2000000000000000      ,    0.0000000000000000       },
    {   -42.899999999999999      ,   -17.800000000000001      ,   0.40000000000000002      ,    0.0000000000000000       },
    {    2.6000000000000001      ,    34.299999999999997      ,   0.80000000000000004      ,    0.0000000000000000       },
    {  -0.80000000000000004      ,   -48.600000000000001      ,    2.3999999999999999      ,  -0.10000000000000001       },
    {   -4.9000000000000004      ,    30.500000000000000      ,    3.7000000000000002      ,   0.69999999999999996       },
    {    0.0000000000000000      ,   -43.600000000000001      ,    2.1000000000000001      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -25.399999999999999      ,    1.2000000000000000      ,    0.0000000000000000       },
    {    2.0000000000000000      ,    40.899999999999999      ,   -2.0000000000000000      ,    0.0000000000000000       },
    {   -2.1000000000000001      ,    26.100000000000001      ,   0.59999999999999998      ,    0.0000000000000000       },
    {    22.600000000000001      ,   -3.2000000000000002      ,  -0.50000000000000000      ,  -0.50000000000000000       },
    {   -7.5999999999999996      ,    24.899999999999999      ,  -0.40000000000000002      ,  -0.20000000000000001       },
    {   -6.2000000000000002      ,    34.899999999999999      ,    1.7000000000000000      ,   0.29999999999999999       },
    {    2.0000000000000000      ,    17.399999999999999      ,  -0.40000000000000002      ,   0.10000000000000001       },
    {   -3.8999999999999999      ,    20.500000000000000      ,    2.3999999999999999      ,   0.59999999999999998       }
  };
  double eps[         194 ][4] = {
    {    9205365.8000000007      ,   -1506.2000000000000      ,    885.70000000000005      ,  -0.20000000000000001       },
    {    573095.90000000002      ,   -570.20000000000005      ,   -305.00000000000000      ,  -0.29999999999999999       },
    {    97845.500000000000      ,    147.80000000000001      ,   -48.799999999999997      ,  -0.20000000000000001       },
    {   -89753.600000000006      ,    28.000000000000000      ,    46.899999999999999      ,    0.0000000000000000       },
    {    7406.6999999999998      ,   -327.10000000000002      ,   -18.199999999999999      ,   0.80000000000000004       },
    {    22442.299999999999      ,   -22.300000000000001      ,   -67.599999999999994      ,    0.0000000000000000       },
    {   -683.60000000000002      ,    46.799999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    20070.700000000001      ,    36.000000000000000      ,    1.6000000000000001      ,    0.0000000000000000       },
    {    12893.799999999999      ,    39.500000000000000      ,   -6.2000000000000002      ,    0.0000000000000000       },
    {   -9593.2000000000007      ,    14.400000000000000      ,    30.199999999999999      ,  -0.10000000000000001       },
    {   -6899.5000000000000      ,    4.7999999999999998      ,  -0.59999999999999998      ,    0.0000000000000000       },
    {   -5332.5000000000000      ,  -0.10000000000000001      ,    2.7000000000000002      ,    0.0000000000000000       },
    {   -125.20000000000000      ,    10.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -3323.4000000000001      ,  -0.90000000000000002      ,  -0.29999999999999999      ,    0.0000000000000000       },
    {    3142.3000000000002      ,    8.9000000000000004      ,   0.29999999999999999      ,    0.0000000000000000       },
    {    2552.5000000000000      ,    7.2999999999999998      ,   -1.2000000000000000      ,    0.0000000000000000       },
    {    2634.4000000000001      ,    8.8000000000000007      ,   0.20000000000000001      ,    0.0000000000000000       },
    {   -2424.4000000000001      ,    1.6000000000000001      ,  -0.40000000000000002      ,    0.0000000000000000       },
    {   -123.30000000000000      ,    3.8999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    1642.4000000000001      ,    7.2999999999999998      ,  -0.80000000000000004      ,    0.0000000000000000       },
    {    47.899999999999999      ,    3.2000000000000002      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    1321.2000000000000      ,    6.2000000000000002      ,  -0.59999999999999998      ,    0.0000000000000000       },
    {   -1234.0999999999999      ,  -0.29999999999999999      ,   0.59999999999999998      ,    0.0000000000000000       },
    {   -1076.5000000000000      ,  -0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -61.600000000000001      ,    1.8000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -55.399999999999999      ,    1.6000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    856.89999999999998      ,   -4.9000000000000004      ,   -2.1000000000000001      ,    0.0000000000000000       },
    {   -800.70000000000005      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    685.10000000000002      ,  -0.59999999999999998      ,   -3.7999999999999998      ,    0.0000000000000000       },
    {   -16.899999999999999      ,   -1.5000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    695.70000000000005      ,    1.8000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    642.20000000000005      ,   -2.6000000000000001      ,   -1.6000000000000001      ,    0.0000000000000000       },
    {    13.300000000000001      ,    1.1000000000000001      ,  -0.10000000000000001      ,    0.0000000000000000       },
    {    521.89999999999998      ,    1.6000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    325.80000000000001      ,    2.0000000000000000      ,  -0.10000000000000001      ,    0.0000000000000000       },
    {   -325.10000000000002      ,  -0.50000000000000000      ,   0.90000000000000002      ,    0.0000000000000000       },
    {    10.100000000000000      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    334.50000000000000      ,    1.6000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    307.10000000000002      ,   0.40000000000000002      ,  -0.90000000000000002      ,    0.0000000000000000       },
    {    327.19999999999999      ,   0.50000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -304.60000000000002      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    304.00000000000000      ,   0.59999999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -276.80000000000001      ,  -0.50000000000000000      ,   0.10000000000000001      ,    0.0000000000000000       },
    {    268.89999999999998      ,    1.3000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    271.80000000000001      ,    1.1000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    271.50000000000000      ,  -0.40000000000000002      ,  -0.80000000000000004      ,    0.0000000000000000       },
    {   -5.2000000000000002      ,   0.50000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -220.50000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -20.100000000000001      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -191.00000000000000      ,   0.10000000000000001      ,   0.50000000000000000      ,    0.0000000000000000       },
    {   -4.0999999999999996      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    130.59999999999999      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    3.0000000000000000      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    122.90000000000001      ,   0.80000000000000004      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    3.7000000000000002      ,  -0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    123.09999999999999      ,   0.40000000000000002      ,  -0.29999999999999999      ,    0.0000000000000000       },
    {   -52.700000000000003      ,    15.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    120.70000000000000      ,   0.29999999999999999      ,  -0.29999999999999999      ,    0.0000000000000000       },
    {    4.0000000000000000      ,  -0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    126.50000000000000      ,   0.50000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    112.70000000000000      ,   0.50000000000000000      ,  -0.29999999999999999      ,    0.0000000000000000       },
    {   -106.09999999999999      ,  -0.29999999999999999      ,   0.29999999999999999      ,    0.0000000000000000       },
    {   -112.90000000000001      ,  -0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    3.6000000000000001      ,  -0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    107.40000000000001      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -10.900000000000000      ,   0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.90000000000000002      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    85.400000000000006      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -88.799999999999997      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -71.000000000000000      ,  -0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -70.299999999999997      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    64.500000000000000      ,   0.40000000000000002      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    69.799999999999997      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    66.099999999999994      ,   0.40000000000000002      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -61.000000000000000      ,  -0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -59.500000000000000      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -55.600000000000001      ,    0.0000000000000000      ,   0.20000000000000001      ,    0.0000000000000000       },
    {    51.700000000000003      ,   0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -49.000000000000000      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -52.700000000000003      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -49.600000000000001      ,    1.3999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    46.299999999999997      ,   0.40000000000000002      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    49.600000000000001      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -5.0999999999999996      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -44.000000000000000      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -39.899999999999999      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -39.500000000000000      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -3.8999999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -42.100000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -17.199999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -2.2999999999999998      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -39.200000000000003      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -38.399999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    36.799999999999997      ,   0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    34.600000000000001      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -32.700000000000003      ,   0.29999999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    30.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.40000000000000002      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    29.300000000000001      ,   0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    31.600000000000001      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.80000000000000004      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -27.899999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    2.8999999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -25.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    25.000000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    27.500000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -24.399999999999999      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    24.899999999999999      ,   0.20000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -22.800000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.90000000000000002      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    24.399999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    23.899999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    22.500000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    20.800000000000001      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    20.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    21.500000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -20.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    1.3999999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.20000000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    19.000000000000000      ,    0.0000000000000000      ,  -0.10000000000000001      ,    0.0000000000000000       },
    {    20.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -2.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -17.600000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    19.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -2.3999999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -18.399999999999999      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    17.100000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.40000000000000002      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    18.399999999999999      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,    17.399999999999999      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.59999999999999998      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -15.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -16.800000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    16.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -2.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.5000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -14.300000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    14.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -13.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -14.300000000000001      ,  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -13.699999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    13.100000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.7000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -12.800000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   -14.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    12.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -12.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {  -0.80000000000000004      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    10.900000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -10.800000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    10.500000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -10.400000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -11.199999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    10.500000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -1.3999999999999999      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    0.0000000000000000      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.69999999999999996      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -10.300000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -10.000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    9.5999999999999996      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {    9.4000000000000004      ,   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.59999999999999998      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -87.700000000000003      ,    4.4000000000000004      ,  -0.40000000000000002      ,   -6.2999999999999998       },
    {    46.299999999999997      ,    22.399999999999999      ,   0.50000000000000000      ,   -2.3999999999999999       },
    {    15.600000000000000      ,   -3.3999999999999999      ,   0.10000000000000001      ,   0.40000000000000002       },
    {    5.2000000000000002      ,    5.7999999999999998      ,   0.20000000000000001      ,  -0.10000000000000001       },
    {   -30.100000000000001      ,    26.899999999999999      ,   0.69999999999999996      ,    0.0000000000000000       },
    {    23.199999999999999      ,  -0.50000000000000000      ,    0.0000000000000000      ,   0.59999999999999998       },
    {    1.0000000000000000      ,    23.199999999999999      ,    3.3999999999999999      ,    0.0000000000000000       },
    {   -12.199999999999999      ,   -4.2999999999999998      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -2.1000000000000001      ,   -3.7000000000000002      ,  -0.20000000000000001      ,   0.10000000000000001       },
    {   -18.600000000000001      ,   -3.7999999999999998      ,  -0.40000000000000002      ,    1.8000000000000000       },
    {    5.5000000000000000      ,   -18.699999999999999      ,   -1.8000000000000000      ,  -0.50000000000000000       },
    {   -5.5000000000000000      ,   -18.699999999999999      ,    1.8000000000000000      ,  -0.50000000000000000       },
    {    18.399999999999999      ,   -3.6000000000000001      ,   0.29999999999999999      ,   0.90000000000000002       },
    {  -0.59999999999999998      ,    1.3000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -5.5999999999999996      ,   -19.500000000000000      ,    1.8999999999999999      ,    0.0000000000000000       },
    {    5.5000000000000000      ,   -19.100000000000001      ,   -1.8999999999999999      ,    0.0000000000000000       },
    {   -17.300000000000001      ,  -0.80000000000000004      ,    0.0000000000000000      ,   0.90000000000000002       },
    {   -3.2000000000000002      ,   -8.3000000000000007      ,  -0.80000000000000004      ,   0.29999999999999999       },
    {  -0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -5.4000000000000004      ,    7.7999999999999998      ,  -0.29999999999999999      ,    0.0000000000000000       },
    {   -14.800000000000001      ,    1.3999999999999999      ,    0.0000000000000000      ,   0.29999999999999999       },
    {   -3.7999999999999998      ,   0.40000000000000002      ,    0.0000000000000000      ,  -0.20000000000000001       },
    {    12.600000000000000      ,    3.2000000000000002      ,   0.50000000000000000      ,   -1.5000000000000000       },
    {   0.10000000000000001      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -13.600000000000000      ,    2.3999999999999999      ,  -0.10000000000000001      ,    0.0000000000000000       },
    {   0.90000000000000002      ,    1.2000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   -11.900000000000000      ,  -0.50000000000000000      ,    0.0000000000000000      ,   0.29999999999999999       },
    {   0.40000000000000002      ,    12.000000000000000      ,   0.29999999999999999      ,  -0.20000000000000001       },
    {    8.3000000000000007      ,    6.0999999999999996      ,  -0.10000000000000001      ,   0.10000000000000001       },
    {    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000      ,    0.0000000000000000       },
    {   0.40000000000000002      ,   -10.800000000000001      ,   0.29999999999999999      ,    0.0000000000000000       },
    {    9.5999999999999996      ,    2.2000000000000002      ,   0.29999999999999999      ,   -1.2000000000000000       }
  };

  /*  interval between fundamental epoch j2000.0 and given epoch (jc). */
  t = (date-DJM0)/DJC;

  /*  mean anomaly of the moon. */
  el  = 134.96340251*PAL__DD2R+
    fmod(t*(1717915923.2178+
            t*(        31.8792+
                       t*(         0.051635+
                                   t*(       - 0.00024470)))),TURNAS)*PAL__DAS2R;

  /*  mean anomaly of the sun. */
  elp = 357.52910918*PAL__DD2R+
    fmod(t*( 129596581.0481+
             t*(       - 0.5532+
                       t*(         0.000136+
                                   t*(       - 0.00001149)))),TURNAS)*PAL__DAS2R;
      
  /*  mean argument of the latitude of the moon. */
  f   =  93.27209062*PAL__DD2R+
    fmod(t*(1739527262.8478+
           t*(      - 12.7512+
                    t*(      -  0.001037+
                             t*(         0.00000417)))),TURNAS)*PAL__DAS2R;

  /*  mean elongation of the moon from the sun. */
  d   = 297.85019547*PAL__DD2R+
    fmod(t*(1602961601.2090+
           t*(       - 6.3706+
                     t*(         0.006539+
                                 t*(       - 0.00003169)))),TURNAS)*PAL__DAS2R;

  /*  mean longitude of the ascending node of the moon. */
  om  = 125.04455501*PAL__DD2R+
    fmod(t*( - 6962890.5431+
            t*(         7.4722+
                        t*(         0.007702+
                                    t*(       - 0.00005939)))),TURNAS)*PAL__DAS2R;

  /*  mean longitude of venus. */
  ve    = 181.97980085*PAL__DD2R+fmod(210664136.433548*t,TURNAS)*PAL__DAS2R;

  /*  mean longitude of mars.*/
  ma    = 355.43299958*PAL__DD2R+fmod( 68905077.493988*t,TURNAS)*PAL__DAS2R;

  /*  mean longitude of jupiter. */
  ju    =  34.35151874*PAL__DD2R+fmod( 10925660.377991*t,TURNAS)*PAL__DAS2R;

  /*  mean longitude of saturn. */
  sa    =  50.07744430*PAL__DD2R+fmod(  4399609.855732*t,TURNAS)*PAL__DAS2R;

  /*  geodesic nutation (fukushima 1991) in microarcsec. */
  dp = -153.1*sin(elp)-1.9*sin(2*elp);
  de = 0.0;

  /*  shirai & fukushima (2001) nutation series. */
  for (j=NTERMS-1; j >= 0; j--) {
    theta = ((double)na[j][0])*el+
      ((double)na[j][1])*elp+
      ((double)na[j][2])*f+
      ((double)na[j][3])*d+
      ((double)na[j][4])*om+
      ((double)na[j][5])*ve+
      ((double)na[j][6])*ma+
      ((double)na[j][7])*ju+
      ((double)na[j][8])*sa;
    c = cos(theta);
    s = sin(theta);
    dp += (psi[j][0] + psi[j][2]*t)*c + (psi[j][1] + psi[j][3]*t)*s;
    de += (eps[j][0] + eps[j][2]*t)*c + (eps[j][1] + eps[j][3]*t)*s;
  }

  /*  change of units, and addition of the precession correction.*/
  *dpsi = (dp*1e-6-0.042888-0.29856*t)*PAL__DAS2R;
  *deps = (de*1e-6-0.005171-0.02408*t)*PAL__DAS2R;

  /*  mean obliquity of date (simon et al. 1994). */
  *eps0 = (84381.412+
          (-46.80927+
           (-0.000152+
            (0.0019989+
             (-0.00000051+
              (-0.000000025)*t)*t)*t)*t)*t)*PAL__DAS2R;

}
