/*
*+
*  Name:
*     palRefro

*  Purpose:
*     Atmospheric refraction for radio and optical/IR wavelengths

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palRefro( double zobs, double hm, double tdk, double pmb,
*                    double rh, double wl, double phi, double tlr,
*                    double eps, double * ref ) {

*  Arguments:
*     zobs = double (Given)
*        Observed zenith distance of the source (radian)
*     hm = double (Given)
*        Height of the observer above sea level (metre)
*     tdk = double (Given)
*        Ambient temperature at the observer (K)
*     pmb = double (Given)
*        Pressure at the observer (millibar)
*     rh = double (Given)
*        Relative humidity at the observer (range 0-1)
*     wl = double (Given)
*        Effective wavelength of the source (micrometre)
*     phi = double (Given)
*        Latitude of the observer (radian, astronomical)
*     tlr = double (Given)
*        Temperature lapse rate in the troposphere (K/metre)
*     eps = double (Given)
*        Precision required to terminate iteration (radian)
*     ref = double * (Returned)
*        Refraction: in vacuao ZD minus observed ZD (radian)

*  Description:
*     Calculates the atmospheric refraction for radio and optical/IR
*     wavelengths.

*  Authors:
*     PTW: Patrick T. Wallace
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - A suggested value for the TLR argument is 0.0065.  The
*     refraction is significantly affected by TLR, and if studies
*     of the local atmosphere have been carried out a better TLR
*     value may be available.  The sign of the supplied TLR value
*     is ignored.
*
*     - A suggested value for the EPS argument is 1E-8.  The result is
*     usually at least two orders of magnitude more computationally
*     precise than the supplied EPS value.
*
*     - The routine computes the refraction for zenith distances up
*     to and a little beyond 90 deg using the method of Hohenkerk
*     and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
*     in the Explanatory Supplement, 1992 edition - see section 3.281).
*
*     - The code is a development of the optical/IR refraction subroutine
*     AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
*     support the radio case.  Apart from merely cosmetic changes, the
*     following modifications to the original HMNAO optical/IR refraction
*     code have been made:
*
*     .  The angle arguments have been changed to radians.
*
*     .  Any value of ZOBS is allowed (see note 6, below).
*
*     .  Other argument values have been limited to safe values.
*
*     .  Murray's values for the gas constants have been used
*        (Vectorial Astrometry, Adam Hilger, 1983).
*
*     .  The numerical integration phase has been rearranged for
*        extra clarity.
*
*     .  A better model for Ps(T) has been adopted (taken from
*        Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
*
*     .  More accurate expressions for Pwo have been adopted
*        (again from Gill 1982).
*
*     .  The formula for the water vapour pressure, given the
*        saturation pressure and the relative humidity, is from
*        Crane (1976), expression 2.5.5.

*     .  Provision for radio wavelengths has been added using
*        expressions devised by A.T.Sinclair, RGO (private
*        communication 1989).  The refractivity model currently
*        used is from J.M.Rueger, "Refractive Index Formulae for
*        Electronic Distance Measurement with Radio and Millimetre
*        Waves", in Unisurv Report S-68 (2002), School of Surveying
*        and Spatial Information Systems, University of New South
*        Wales, Sydney, Australia.
*
*     .  The optical refractivity for dry air is from Resolution 3 of
*        the International Association of Geodesy adopted at the XXIIth
*        General Assembly in Birmingham, UK, 1999.
*
*     .  Various small changes have been made to gain speed.
*
*     - The radio refraction is chosen by specifying WL > 100 micrometres.
*     Because the algorithm takes no account of the ionosphere, the
*     accuracy deteriorates at low frequencies, below about 30 MHz.
*
*     - Before use, the value of ZOBS is expressed in the range +/- pi.
*     If this ranged ZOBS is -ve, the result REF is computed from its
*     absolute value before being made -ve to match.  In addition, if
*     it has an absolute value greater than 93 deg, a fixed REF value
*     equal to the result for ZOBS = 93 deg is returned, appropriately
*     signed.
*
*     - As in the original Hohenkerk and Sinclair algorithm, fixed values
*     of the water vapour polytrope exponent, the height of the
*     tropopause, and the height at which refraction is negligible are
*     used.
*
*     - The radio refraction has been tested against work done by
*     Iain Coulson, JACH, (private communication 1995) for the
*     James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
*     agreement at the 0.1 arcsec level is achieved for moderate ZD,
*     worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
*     humid sea-level sites the accuracy will not be as good.
*
*     - It should be noted that the relative humidity RH is formally
*     defined in terms of "mixing ratio" rather than pressures or
*     densities as is often stated.  It is the mass of water per unit
*     mass of dry air divided by that for saturated air at the same
*     temperature and pressure (see Gill 1982).

*     - The algorithm is designed for observers in the troposphere. The
*     supplied temperature, pressure and lapse rate are assumed to be
*     for a point in the troposphere and are used to define a model
*     atmosphere with the tropopause at 11km altitude and a constant
*     temperature above that.  However, in practice, the refraction
*     values returned for stratospheric observers, at altitudes up to
*     25km, are quite usable.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version, direct port of SLA Fortran source.
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

#include <math.h>

#include "pal.h"
#include "pal1.h"
#include "palmac.h"

void palRefro( double zobs, double hm, double tdk, double pmb,
               double rh, double wl, double phi, double tlr,
               double eps, double * ref ) {

  /*
   *  Fixed parameters
   */

  /*  93 degrees in radians */
  const double D93 = 1.623156204;
  /*  Universal gas constant */
  const double GCR = 8314.32;
  /*  Molecular weight of dry air */
  const double DMD = 28.9644;
  /*  Molecular weight of water vapour */
  const double DMW = 18.0152;
  /*  Mean Earth radius (metre) */
  const double S = 6378120.;
  /*  Exponent of temperature dependence of water vapour pressure */
  const double DELTA = 18.36;
  /*  Height of tropopause (metre) */
  const double HT = 11000.;
  /*  Upper limit for refractive effects (metre) */
  const double HS = 80000.;
  /*  Numerical integration: maximum number of strips. */
  const int ISMAX=16384l;

  /* Local variables */
  int is, k, n, i, j;

  int optic, loop; /* booleans */

  double zobs1,zobs2,hmok,tdkok,pmbok,rhok,wlok,alpha,
    tol,wlsq,gb,a,gamal,gamma,gamm2,delm2,
    tdc,psat,pwo,w,
    c1,c2,c3,c4,c5,c6,r0,tempo,dn0,rdndr0,sk0,f0,
    rt,tt,dnt,rdndrt,sine,zt,ft,dnts,rdndrp,zts,fts,
    rs,dns,rdndrs,zs,fs,refold,z0,zrange,fb,ff,fo,fe,
    h,r,sz,rg,dr,tg,dn,rdndr,t,f,refp,reft;

  /*  The refraction integrand */
#define refi(DN,RDNDR) RDNDR/(DN+RDNDR)

  /*  Transform ZOBS into the normal range. */
  zobs1 = palDrange(zobs);
  zobs2 = DMIN(fabs(zobs1),D93);

  /*  keep other arguments within safe bounds. */
  hmok = DMIN(DMAX(hm,-1e3),HS);
  tdkok = DMIN(DMAX(tdk,100.0),500.0);
  pmbok = DMIN(DMAX(pmb,0.0),10000.0);
  rhok = DMIN(DMAX(rh,0.0),1.0);
  wlok = DMAX(wl,0.1);
  alpha = DMIN(DMAX(fabs(tlr),0.001),0.01);

  /*  tolerance for iteration. */
  tol = DMIN(DMAX(fabs(eps),1e-12),0.1)/2.0;

  /*  decide whether optical/ir or radio case - switch at 100 microns. */
  optic = wlok < 100.0;

  /*  set up model atmosphere parameters defined at the observer. */
  wlsq = wlok*wlok;
  gb = 9.784*(1.0-0.0026*cos(phi+phi)-0.00000028*hmok);
  if (optic) {
    a = (287.6155+(1.62887+0.01360/wlsq)/wlsq) * 273.15e-6/1013.25;
  } else {
    a = 77.6890e-6;
  }
  gamal = (gb*DMD)/GCR;
  gamma = gamal/alpha;
  gamm2 = gamma-2.0;
  delm2 = DELTA-2.0;
  tdc = tdkok-273.15;
  psat = pow(10.0,(0.7859+0.03477*tdc)/(1.0+0.00412*tdc)) *
    (1.0+pmbok*(4.5e-6+6.0e-10*tdc*tdc));
  if (pmbok > 0.0) {
    pwo = rhok*psat/(1.0-(1.0-rhok)*psat/pmbok);
  } else {
    pwo = 0.0;
  }
  w = pwo*(1.0-DMW/DMD)*gamma/(DELTA-gamma);
  c1 = a*(pmbok+w)/tdkok;
  if (optic) {
    c2 = (a*w+11.2684e-6*pwo)/tdkok;
  } else {
    c2 = (a*w+6.3938e-6*pwo)/tdkok;
  }
  c3 = (gamma-1.0)*alpha*c1/tdkok;
  c4 = (DELTA-1.0)*alpha*c2/tdkok;
  if (optic) {
    c5 = 0.0;
    c6 = 0.0;
  } else {
    c5 = 375463e-6*pwo/tdkok;
    c6 = c5*delm2*alpha/(tdkok*tdkok);
  }

  /*  conditions at the observer. */
  r0 = S+hmok;
  pal1Atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,
           r0,&tempo,&dn0,&rdndr0);
  sk0 = dn0*r0*sin(zobs2);
  f0 = refi(dn0,rdndr0);

  /*  conditions in the troposphere at the tropopause. */
  rt = S+DMAX(HT,hmok);
  pal1Atmt(r0,tdkok,alpha,gamm2,delm2,c1,c2,c3,c4,c5,c6,
           rt,&tt,&dnt,&rdndrt);
  sine = sk0/(rt*dnt);
  zt = atan2(sine,sqrt(DMAX(1.0-sine*sine,0.0)));
  ft = refi(dnt,rdndrt);

  /*  conditions in the stratosphere at the tropopause. */
  pal1Atms(rt,tt,dnt,gamal,rt,&dnts,&rdndrp);
  sine = sk0/(rt*dnts);
  zts = atan2(sine,sqrt(DMAX(1.0-sine*sine,0.0)));
  fts = refi(dnts,rdndrp);

  /*  conditions at the stratosphere limit. */
  rs = S+HS;
  pal1Atms(rt,tt,dnt,gamal,rs,&dns,&rdndrs);
  sine = sk0/(rs*dns);
  zs = atan2(sine,sqrt(DMAX(1.0-sine*sine,0.0)));
  fs = refi(dns,rdndrs);

  /*  variable initialization to avoid compiler warning. */
  reft = 0.0;

  /*  integrate the refraction integral in two parts;  first in the
   *  troposphere (k=1), then in the stratosphere (k=2). */

  for (k=1; k<=2; k++) {

    /*     initialize previous refraction to ensure at least two iterations. */
    refold = 1.0;

    /*     start off with 8 strips. */
    is = 8;

    /*     start z, z range, and start and end values. */
    if (k==1) {
      z0 = zobs2;
      zrange = zt-z0;
      fb = f0;
      ff = ft;
    } else {
      z0 = zts;
      zrange = zs-z0;
      fb = fts;
      ff = fs;
    }

    /*     sums of odd and even values. */
    fo = 0.0;
    fe = 0.0;

    /*     first time through the loop we have to do every point. */
    n = 1;

    /*     start of iteration loop (terminates at specified precision). */
    loop = 1;
    while (loop) {

      /*        strip width. */
      h = zrange/((double)is);

      /*        initialize distance from earth centre for quadrature pass. */
      if (k == 1) {
        r = r0;
      } else {
        r = rt;
      }

      /*        one pass (no need to compute evens after first time). */
      for (i=1; i<is; i+=n) {

              /*           sine of observed zenith distance. */
        sz = sin(z0+h*(double)(i));

        /*           find r (to the nearest metre, maximum four iterations). */
        if (sz > 1e-20) {
          w = sk0/sz;
          rg = r;
          dr = 1.0e6;
          j = 0;
          while ( fabs(dr) > 1.0 && j < 4 ) {
            j++;
            if (k==1) {
              pal1Atmt(r0,tdkok,alpha,gamm2,delm2,
                       c1,c2,c3,c4,c5,c6,rg,&tg,&dn,&rdndr);
            } else {
              pal1Atms(rt,tt,dnt,gamal,rg,&dn,&rdndr);
            }
            dr = (rg*dn-w)/(dn+rdndr);
            rg = rg-dr;
          }
          r = rg;
        }

        /*           find the refractive index and integrand at r. */
        if (k==1) {
          pal1Atmt(r0,tdkok,alpha,gamm2,delm2,
                   c1,c2,c3,c4,c5,c6,r,&t,&dn,&rdndr);
        } else {
          pal1Atms(rt,tt,dnt,gamal,r,&dn,&rdndr);
        }
        f = refi(dn,rdndr);

        /*           accumulate odd and (first time only) even values. */
        if (n==1 && i%2 == 0) {
          fe += f;
        } else {
          fo += f;
        }
      }

      /*        evaluate the integrand using simpson's rule. */
      refp = h*(fb+4.0*fo+2.0*fe+ff)/3.0;

      /*        has the required precision been achieved (or can't be)? */
      if (fabs(refp-refold) > tol && is < ISMAX) {

        /*           no: prepare for next iteration.*/

        /*           save current value for convergence test. */
        refold = refp;

        /*           double the number of strips. */
        is += is;

        /*           sum of all current values = sum of next pass's even values. */
        fe += fo;

        /*           prepare for new odd values. */
        fo = 0.0;

        /*           skip even values next time. */
        n = 2;
      } else {

        /*           yes: save troposphere component and terminate the loop. */
        if (k==1) reft = refp;
        loop = 0;
      }
    }
  }

  /*  result. */
  *ref = reft+refp;
  if (zobs1 < 0.0) *ref = -(*ref);

}
