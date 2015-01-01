/*
*+
*  Name:
*     palUnpcd

*  Purpose:
*     Remove pincushion/barrel distortion

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palUnpcd( double disco, double * x, double * y );

*  Arguments:
*     disco = double (Given)
*        Pincushion/barrel distortion coefficient.
*     x = double * (Given & Returned)
*        On input the distorted X coordinate, on output
*        the tangent-plane X coordinate.
*     y = double * (Given & Returned)
*        On input the distorted Y coordinate, on output
*        the tangent-plane Y coordinate.

*  Description:
*     Remove pincushion/barrel distortion from a distorted [x,y] to give
*     tangent-plane [x,y].

*  Authors:
*     PTW: Pat Wallace (RAL)
*     TIMJ: Tim Jenness (Cornell)
*     {enter_new_authors_here}

*  Notes:
*     - The distortion is of the form RP = R*(1+C*R^2), where R is
*       the radial distance from the tangent point, C is the DISCO
*       argument, and RP is the radial distance in the presence of
*       the distortion.
*
*     - For pincushion distortion, C is +ve;  for barrel distortion,
*       C is -ve.
*
*     - For X,Y in "radians" - units of one projection radius,
*       which in the case of a photograph is the focal length of
*       the camera - the following DISCO values apply:
*
*           Geometry          DISCO
*
*           astrograph         0.0
*           Schmidt           -0.3333
*           AAT PF doublet  +147.069
*           AAT PF triplet  +178.585
*           AAT f/8          +21.20
*           JKT f/8          +13.32
*
*     - The present routine is a rigorous inverse of the companion
*       routine palPcd.  The expression for RP in Note 1 is rewritten
*       in the form x^3+a*x+b=0 and solved by standard techniques.
*
*     - Cases where the cubic has multiple real roots can sometimes
*       occur, corresponding to extreme instances of barrel distortion
*       where up to three different undistorted [X,Y]s all produce the
*       same distorted [X,Y].  However, only one solution is returned,
*       the one that produces the smallest change in [X,Y].

*  See Also:
*     palPcd

*  History:
*     2000-09-03 (PTW):
*        SLALIB implementation.
*     2015-01-01 (TIMJ):
*        Initial version
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2000 Rutherford Appleton Laboratory.
*     Copyright (C) 2015 Cornell University
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
*     along with this program.  If not, see <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "pal.h"
#include "palmac.h"

/* copysign is C99 */
#if HAVE_COPYSIGN
# define COPYSIGN copysign
#else
# define COPYSIGN(a,b) DSIGN(a,b)
#endif

void palUnpcd( double disco, double * x, double *y ) {

  const double THIRD = 1.0/3.0;

  double rp,q,r,d,w,s,t,f,c,t3,f1,f2,f3,w1,w2,w3;
  double c2;

  /*  Distance of the point from the origin. */
  rp = sqrt( (*x)*(*x)+(*y)*(*y));

  /*  If zero, or if no distortion, no action is necessary. */
  if (rp != 0.0 && disco != 0.0) {

    /*     Begin algebraic solution. */
    q = 1.0/(3.0*disco);
    r = rp/(2.0*disco);
    w = q*q*q+r*r;

    /* Continue if one real root, or three of which only one is positive. */
    if (w > 0.0) {

      d = sqrt(w);
      w = r+d;
      s = COPYSIGN(pow(fabs(w),THIRD),w);
      w = r-d;
      t = COPYSIGN(pow(fabs(w),THIRD),w);
      f = s+t;

    } else {
      /* Three different real roots:  use geometrical method instead. */
      w = 2.0/sqrt(-3.0*disco);
      c = 4.0*rp/(disco*w*w*w);
      c2 = c*c;
      s = sqrt(1.0-DMIN(c2,1.0));
      t3 = atan2(s,c);

      /* The three solutions. */
      f1 = w*cos((PAL__D2PI-t3)/3.0);
      f2 = w*cos((t3)/3.0);
      f3 = w*cos((PAL__D2PI+t3)/3.0);

      /* Pick the one that moves [X,Y] least. */
      w1 = fabs(f1-rp);
      w2 = fabs(f2-rp);
      w3 = fabs(f3-rp);
      if (w1 < w2) {
        f = ( w1 < w3 ? f1 : f3 );
      } else {
        f = ( w2 < w3 ? f2 : f3 );
      }
    }

    /* Remove the distortion. */
    f = f/rp;
    *x *= f;
    *y *= f;
  }
}
