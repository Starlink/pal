/*
*+
*  Name:
*     palAtmdsp

*  Purpose:
*     Apply atmospheric-dispersion adjustments to refraction coefficients

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palAtmdsp( double tdk, double pmb, double rh, double wl1,
*                     double a1, double b1, double wl2, double *a2, double *b2 );


*  Arguments:
*     tdk = double (Given)
*        Ambient temperature, K
*     pmb = double (Given)
*        Ambient pressure, millibars
*     rh = double (Given)
*        Ambient relative humidity, 0-1
*     wl1 = double (Given)
*        Reference wavelength, micrometre (0.4 recommended)
*     a1 = double (Given)
*        Refraction coefficient A for wavelength wl1 (radians)
*     b1 = double (Given)
*        Refraction coefficient B for wavelength wl1 (radians)
*     wl2 = double (Given)
*        Wavelength for which adjusted A,B required
*     a2 = double * (Returned)
*        Refraction coefficient A for wavelength WL2 (radians)
*     b2 = double * (Returned)
*        Refraction coefficient B for wavelength WL2 (radians)

*  Description:
*     Apply atmospheric-dispersion adjustments to refraction coefficients.

*  Authors:
*     TIMJ: Tim Jenness
*     PTW: Patrick Wallace
*     {enter_new_authors_here}

*  Notes:
*     - To use this routine, first call palRefco specifying WL1 as the
*     wavelength.  This yields refraction coefficients A1,B1, correct
*     for that wavelength.  Subsequently, calls to palAtmdsp specifying
*     different wavelengths will produce new, slightly adjusted
*     refraction coefficients which apply to the specified wavelength.
*
*     - Most of the atmospheric dispersion happens between 0.7 micrometre
*     and the UV atmospheric cutoff, and the effect increases strongly
*     towards the UV end.  For this reason a blue reference wavelength
*     is recommended, for example 0.4 micrometres.
*
*     - The accuracy, for this set of conditions:
*
*        height above sea level    2000 m
*                      latitude    29 deg
*                      pressure    793 mb
*                   temperature    17 degC
*                      humidity    50%
*                    lapse rate    0.0065 degC/m
*          reference wavelength    0.4 micrometre
*                star elevation    15 deg
*
*     is about 2.5 mas RMS between 0.3 and 1.0 micrometres, and stays
*     within 4 mas for the whole range longward of 0.3 micrometres
*     (compared with a total dispersion from 0.3 to 20.0 micrometres
*     of about 11 arcsec).  These errors are typical for ordinary
*     conditions and the given elevation;  in extreme conditions values
*     a few times this size may occur, while at higher elevations the
*     errors become much smaller.
*
*     - If either wavelength exceeds 100 micrometres, the radio case
*     is assumed and the returned refraction coefficients are the
*     same as the given ones.  Note that radio refraction coefficients
*     cannot be turned into optical values using this routine, nor
*     vice versa.
*
*     - The algorithm consists of calculation of the refractivity of the
*     air at the observer for the two wavelengths, using the methods
*     of the palRefro routine, and then scaling of the two refraction
*     coefficients according to classical refraction theory.  This
*     amounts to scaling the A coefficient in proportion to (n-1) and
*     the B coefficient almost in the same ratio (see R.M.Green,
*     "Spherical Astronomy", Cambridge University Press, 1985).

*  History:
*     2014-07-15 (TIMJ):
*        Initial version. A direct copy of the Fortran SLA implementation.
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2014 Tim Jenness
*     Copyright (C) 2005 Patrick Wallace
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
#include <math.h>

void palAtmdsp ( double tdk, double pmb, double rh, double wl1,
                 double a1, double b1, double wl2, double *a2, double *b2 ) {

  double f,tdkok,pmbok,rhok;
  double psat,pwo,w1,wlok,wlsq,w2,dn1,dn2;

  /*  Check for radio wavelengths */
  if (wl1 > 100.0 || wl2 > 100.0) {

    /*     Radio: no dispersion */
    *a2 = a1;
    *b2 = b1;

  } else {

    /*     Optical: keep arguments within safe bounds */
    tdkok = DMIN(DMAX(tdk,100.0),500.0);
    pmbok = DMIN(DMAX(pmb,0.0),10000.0);
    rhok = DMIN(DMAX(rh,0.0),1.0);

    /*     Atmosphere parameters at the observer */
    psat = pow(10.0, -8.7115+0.03477*tdkok);
    pwo = rhok*psat;
    w1 = 11.2684e-6*pwo;

    /*     Refractivity at the observer for first wavelength */
    wlok = DMAX(wl1,0.1);
    wlsq = wlok*wlok;
    w2 = 77.5317e-6+(0.43909e-6+0.00367e-6/wlsq)/wlsq;
    dn1 = (w2*pmbok-w1)/tdkok;

    /*     Refractivity at the observer for second wavelength */
    wlok = DMAX(wl2,0.1);
    wlsq = wlok*wlok;
    w2 = 77.5317e-6+(0.43909e-6+0.00367e-6/wlsq)/wlsq;
    dn2 = (w2*pmbok-w1)/tdkok;

    /*     Scale the refraction coefficients (see Green 4.31, p93) */
    if (dn1 != 0.0) {
      f = dn2/dn1;
      *a2 = a1*f;
      *b2 = b1*f;
      if (dn1 != a1) {
        *b2 *= (1.0+dn1*(dn1-dn2)/(2.0*(dn1-a1)));
      }
    }  else {
      *a2 = a1;
      *b2 = b1;
    }
  }

}
