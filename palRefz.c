/*
*+
*  Name:
*     palRefz

*  Purpose:
*     Adjust unrefracted zenith distance

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palRefz ( double zu, double refa, double refb, double *zr );

*  Arguments:
*     zu = double (Given)
*         Unrefracted zenith distance of the source (radians)
*     refa = double (Given)
*         tan Z coefficient (radians)
*     refb = double (Given)
*         tan**3 Z coefficient (radian)
*     zr = double * (Returned)
*         Refracted zenith distance (radians)

*  Description:
*     Adjust an unrefracted zenith distance to include the effect of
*     atmospheric refraction, using the simple A tan Z + B tan**3 Z
*     model (plus special handling for large ZDs).

*  Authors:
*     PTW: Patrick T. Wallace
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - This routine applies the adjustment for refraction in the
*     opposite sense to the usual one - it takes an unrefracted
*     (in vacuo) position and produces an observed (refracted)
*     position, whereas the A tan Z + B tan**3 Z model strictly
*     applies to the case where an observed position is to have the
*     refraction removed.  The unrefracted to refracted case is
*     harder, and requires an inverted form of the text-book
*     refraction models;  the formula used here is based on the
*     Newton-Raphson method.  For the utmost numerical consistency
*     with the refracted to unrefracted model, two iterations are
*     carried out, achieving agreement at the 1D-11 arcseconds level
*     for a ZD of 80 degrees.  The inherent accuracy of the model
*     is, of course, far worse than this - see the documentation for
*     palRefco for more information.
*
*     - At ZD 83 degrees, the rapidly-worsening A tan Z + B tan^3 Z
*     model is abandoned and an empirical formula takes over.  For
*     optical/IR wavelengths, over a wide range of observer heights and
*     corresponding temperatures and pressures, the following levels of
*     accuracy (arcsec, worst case) are achieved, relative to numerical
*     integration through a model atmosphere:
*
*              ZR    error
*
*              80      0.7
*              81      1.3
*              82      2.4
*              83      4.7
*              84      6.2
*              85      6.4
*              86      8
*              87     10
*              88     15
*              89     30
*              90     60
*              91    150         } relevant only to
*              92    400         } high-elevation sites
*
*     For radio wavelengths the errors are typically 50% larger than
*     the optical figures and by ZD 85 deg are twice as bad, worsening
*     rapidly below that.  To maintain 1 arcsec accuracy down to ZD=85
*     at the Green Bank site, Condon (2004) has suggested amplifying
*     the amount of refraction predicted by palRefz below 10.8 deg
*     elevation by the factor (1+0.00195*(10.8-E_t)), where E_t is the
*     unrefracted elevation in degrees.
*
*     The high-ZD model is scaled to match the normal model at the
*     transition point;  there is no glitch.
*
*     - Beyond 93 deg zenith distance, the refraction is held at its
*     93 deg value.
*
*     - See also the routine palRefv, which performs the adjustment in
*     Cartesian Az/El coordinates, and with the emphasis on speed
*     rather than numerical accuracy.

*  References:
*     Condon,J.J., Refraction Corrections for the GBT, PTCS/PN/35.2,
*     NRAO Green Bank, 2004.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version, ported directly from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2004 Rutherford Appleton Laboratory
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

void palRefz ( double zu, double refa, double refb, double *zr ) {

  /* Constants */

  /* Largest usable ZD (deg) */
  const double D93 = 93.0;

  /* ZD at which one model hands over to the other (radians) */
  const double Z83 = 83.0 * PAL__DD2R;

  /* coefficients for high ZD model (used beyond ZD 83 deg) */
  const double C1 = +0.55445;
  const double C2 = -0.01133;
  const double C3 = +0.00202;
  const double C4 = +0.28385;
  const double C5 = +0.02390;

  /* High-ZD-model prefiction (deg) for that point */
  const double REF83 = (C1+C2*7.0+C3*49.0)/(1.0+C4*7.0+C5*49.0);

  double zu1,zl,s,c,t,tsq,tcu,ref,e,e2;

  /*  perform calculations for zu or 83 deg, whichever is smaller */
  zu1 = DMIN(zu,Z83);

  /*  functions of ZD */
  zl = zu1;
  s = sin(zl);
  c = cos(zl);
  t = s/c;
  tsq = t*t;
  tcu = t*tsq;

  /*  refracted zd (mathematically to better than 1 mas at 70 deg) */
  zl = zl-(refa*t+refb*tcu)/(1.0+(refa+3.0*refb*tsq)/(c*c));

  /*  further iteration */
  s = sin(zl);
  c = cos(zl);
  t = s/c;
  tsq = t*t;
  tcu = t*tsq;
  ref = zu1-zl+
    (zl-zu1+refa*t+refb*tcu)/(1.0+(refa+3.0*refb*tsq)/(c*c));

  /*  special handling for large zu */
  if (zu > zu1) {
    e = 90.0-DMIN(D93,zu*PAL__DR2D);
    e2 = e*e;
    ref = (ref/REF83)*(C1+C2*e+C3*e2)/(1.0+C4*e+C5*e2);
  }

  /*  return refracted zd */
  *zr = zu-ref;

}
