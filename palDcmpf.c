/*
*+
*  Name:
*     palDcmpf

*  Purpose:
*     Decompose an [x,y] linear fit into its constituent parameters:
*     zero points, scales, nonperpendicularity and orientation.

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDcmpf ( double coeffs[6], double *xz, double *yz, double *xs,
*                double *ys, double *perp, double *orient )

*  Arguments:
*     coeffs = double[6] (Given)
*        transformation coefficients (see note)
*     xz = double (Returned)
*        x zero point
*     yz = double (Returned)
*        y zero point
*     xs = double (Returned)
*        x scale
*     ys = double (Returned)
*        y scale
*     perp = double (Returned)
*        nonperpendicularity (radians)
*     orient = double (Returned)
*        orientation (radians)

*  Description:
*     The model relates two sets of [x,y] coordinates as follows.
*     Naming the elements of coeffs:
*     ---
*        coeffs[0] = A
*        coeffs[1] = B
*        coeffs[2] = C
*        coeffs[3] = D
*        coeffs[4] = E
*        coeffs[5] = F
*     ---
*     the model transforms coordinates [x1,y1] into coordinates
*     [x2,y2] as follows:
*     ---
*        x2 = A + B * x1 + C * y1
*        y2 = D + E * x1 + F * y1
*     ---
*     The transformation can be decomposed into four steps:
*     ---
*     1) Zero points:
*        x' = xz + x1
*        y' = yz + y1
*
*     2) Scales:
*        x'' = xs * x'
*        y'' = ys * y'
*
*     3) Nonperpendicularity:
*        x''' = cos(perp / 2) * x'' + sin(perp / 2) * y''
*        y''' = sin(perp / 2) * x'' + cos(perp / 2) * y''
*
*     4) Orientation:
*        x2 =  cos(orient) * x''' + sin(orient) * y'''
*        y2 = -sin(orient) * y''' + cos(orient) * y'''
*     ---

*  See also:
*     palFitxy, palPxy, palInvf and palXy2xy

*  Authors:
*     PTW: Pat Wallace (STFC)
*     GSB: Graham Bell (EAO)

*  History:
*     2001-12-19 (PTW):
*        SLALIB implementation.
*     2018-09-13 (GSB):
*        Initial version in C.

*  Copyright:
*     Copyright (C) 2001 Rutherford Appleton Laboratory.
*     Copyright (C) 2018 East Asian Observatory.

*  Licence:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the
*    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
*    Boston, MA  02110-1301  USA
*-
*/

#include <math.h>

#include "pal.h"
#include "palmac.h"

void palDcmpf ( double coeffs[6], double *xz, double *yz, double *xs,
                double *ys, double *perp, double *orient ) {

  double a, b, c, d, e, f, rb2e2, rc2f2, xsc, ysc, p1, p2, p, ws, wc,
    or, hp, shp, chp, sor, cor, det, x0, y0;

  /* Copy the six coefficients. */
  a = coeffs[0];
  b = coeffs[1];
  c = coeffs[2];
  d = coeffs[3];
  e = coeffs[4];
  f = coeffs[5];

  /* Scales. */
  rb2e2 = sqrt(b*b + e*e);
  rc2f2 = sqrt(c*c + f*f);
  if (b*f - c*e >= 0.0) {
    xsc = rb2e2;
  } else {
    b = -b;
    e = -e;
    xsc = -rb2e2;
  }
  ysc = rc2f2;

  /* Non-perpendicularity. */
  if (c != 0.0 || f != 0.0) {
    p1 = atan2(c, f);
  } else {
    p1 = 0.0;
  }
  if (e != 0.0 || b != 0.0) {
    p2 = atan2(e, b);
  } else {
    p2 = 0.0;
  }
  p = palDrange(p1 + p2);

  /* Orientation. */
  ws = c*rb2e2 - e*rc2f2;
  wc = b*rc2f2 + f*rb2e2;
  if (ws != 0.0 || wc != 0.0) {
    or = atan2(ws, wc);
  } else {
    or = 0.0;
  }

  /* Zero points. */
  hp = p / 2.0;
  shp = sin(hp);
  chp = cos(hp);
  sor = sin(or);
  cor = cos(or);
  det = xsc * ysc * (chp + shp) * (chp - shp);
  if (fabs(det) > 0.0) {
    x0 = ysc * (a * (chp*cor - shp*sor) - d * (chp*sor + shp*cor)) / det;
    y0 = xsc * (a * (chp*sor - shp*cor) + d * (chp*cor + shp*sor)) / det;
  } else {
    x0 = 0.0;
    y0 = 0.0;
  }

  /* Results. */
  *xz = x0;
  *yz = y0;
  *xs = xsc;
  *ys = ysc;
  *perp = p;
  *orient = or;
}
