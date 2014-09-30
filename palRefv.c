/*
*+
*  Name:
*     palRefv

*  Purpose:
*     Adjust an unrefracted Cartesian vector to include the effect of atmospheric refraction

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palRefv ( double vu[3], double refa, double refb, double vr[3] );

*  Arguments:
*     vu[3] = double (Given)
*        Unrefracted position of the source (Az/El 3-vector)
*     refa = double (Given)
*        tan Z coefficient (radian)
*     refb = double (Given)
*        tan**3 Z coefficient (radian)
*     vr[3] = double (Returned)
*        Refracted position of the source (Az/El 3-vector)

*  Description:
*     Adjust an unrefracted Cartesian vector to include the effect of
*     atmospheric refraction, using the simple A tan Z + B tan**3 Z
*     model.

*  Authors:
*     TIMJ: Tim Jenness
*     PTW: Patrick Wallace
*     {enter_new_authors_here}

*  Notes:
*     - This routine applies the adjustment for refraction in the
*     opposite sense to the usual one - it takes an unrefracted
*     (in vacuo) position and produces an observed (refracted)
*     position, whereas the A tan Z + B tan**3 Z model strictly
*     applies to the case where an observed position is to have the
*     refraction removed.  The unrefracted to refracted case is
*     harder, and requires an inverted form of the text-book
*     refraction models;  the algorithm used here is equivalent to
*     one iteration of the Newton-Raphson method applied to the above
*     formula.
*
*     - Though optimized for speed rather than precision, the present
*     routine achieves consistency with the refracted-to-unrefracted
*     A tan Z + B tan**3 Z model at better than 1 microarcsecond within
*     30 degrees of the zenith and remains within 1 milliarcsecond to
*     beyond ZD 70 degrees.  The inherent accuracy of the model is, of
*     course, far worse than this - see the documentation for palRefco
*     for more information.
*
*     - At low elevations (below about 3 degrees) the refraction
*     correction is held back to prevent arithmetic problems and
*     wildly wrong results.  For optical/IR wavelengths, over a wide
*     range of observer heights and corresponding temperatures and
*     pressures, the following levels of accuracy (arcsec, worst case)
*     are achieved, relative to numerical integration through a model
*     atmosphere:
*
*              ZD    error
*
*              80      0.7
*              81      1.3
*              82      2.5
*              83      5
*              84     10
*              85     20
*              86     55
*              87    160
*              88    360
*              89    640
*              90   1100
*              91   1700         } relevant only to
*              92   2600         } high-elevation sites
*
*     The results for radio are slightly worse over most of the range,
*     becoming significantly worse below ZD=88 and unusable beyond
*     ZD=90.
*
*     - See also the routine palRefz, which performs the adjustment to
*     the zenith distance rather than in Cartesian Az/El coordinates.
*     The present routine is faster than palRefz and, except very low down,
*     is equally accurate for all practical purposes.  However, beyond
*     about ZD 84 degrees palRefz should be used, and for the utmost
*     accuracy iterative use of palRefro should be considered.

*  History:
*     2014-07-15 (TIMJ):
*        Initial version. A direct copy of the Fortran SLA implementation.
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2014 Tim Jenness
*     Copyright (C) 2004 Patrick Wallace
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

void palRefv ( double vu[3], double refa, double refb, double vr[3] ) {

  double x,y,z1,z,zsq,rsq,r,wb,wt,d,cd,f;

  /*  Initial estimate = unrefracted vector */
  x = vu[0];
  y = vu[1];
  z1 = vu[2];

  /*  Keep correction approximately constant below about 3 deg elevation */
  z = DMAX(z1,0.05);

  /*  One Newton-Raphson iteration */
  zsq = z*z;
  rsq = x*x+y*y;
  r = sqrt(rsq);
  wb = refb*rsq/zsq;
  wt = (refa+wb)/(1.0+(refa+3.0*wb)*(zsq+rsq)/zsq);
  d = wt*r/z;
  cd = 1.0-d*d/2.0;
  f = cd*(1.0-wt);

  /*  Post-refraction x,y,z */
  vr[0] = x*f;
  vr[1] = y*f;
  vr[2] = cd*(z+d*r)+(z1-z);
}
