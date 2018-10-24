/*
*+
*  Name:
*     palPxy

*  Purpose:
*     Given arrays of "expected" and "measured" [x,y] coordinates, and a
*     linear model relating them (as produced by palFitxy), compute
*     the array of "predicted" coordinates and the RMS residuals.

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palPxy ( int np, double xye[][2], double xym[][2], double coeffs[6],
*              double xyp[][2], double *xrms, double *yrms, double *rrms )

*  Arguments:
*     np = int (Given)
*        number of samples
*     xye = double[np][2] (Given)
*        expected [x,y] for each sample
*     xym = double[np][2] (Given)
*        measured [x,y] for each sample
*     coeffs = double[6]
*        coefficients of model (see below)
*     xyp = double[np][2] (Returned)
*        predicted [x,y] for each sample
*     xrms = double (Returned)
*        RMS in x
*     yrms = double (Returned)
*        RMS in y
*     rrms = double (Returned)
*        total RMS (vector sum of xrms and yrms)

*  Description:
*     The model is supplied in the array coeffs.  Naming the
*     elements of coeffs as follows:
*     ---
*        coeffs[0] = A
*        coeffs[1] = B
*        coeffs[2] = C
*        coeffs[3] = D
*        coeffs[4] = E
*        coeffs[5] = F
*     ---
*     the model is applied thus:
*     ---
*        xp = A + B * xm + C * ym
*        yp = D + E * xm + F * ym
*     ---
*     The residuals are (xp - xe) and (yp - ye).
*
*     If np is less than or equal to zero, no coordinates are
*     transformed, and the RMS residuals are all zero.

*  See also:
*     palFitxy, palInvf, palXy2xy and palDcmpf

*  Authors:
*     PTW: Pat Wallace (STFC)
*     GSB: Graham Bell (EAO)

*  History:
*     1996-05-22 (PTW):
*        SLALIB implementation.
*     2018-10-23 (GSB):
*        Initial version in C.

*  Copyright:
*     Copyright (C) 1996 Rutherford Appleton Laboratory
*     Copyright (C) 2018 East Asian Observatory.

*  Licence:
*     This program is free software; you can redistribute it and/or modify
*     it under the terms of the GNU General Public License as published by
*     the Free Software Foundation; either version 2 of the License, or
*     (at your option) any later version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU General Public License for more details.
*
*     You should have received a copy of the GNU General Public License
*     along with this program (see SLA_CONDITIONS); if not, write to the
*     Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
*     Boston, MA  02110-1301  USA
*-
*/

#include <math.h>

#include "pal.h"
#include "palmac.h"

void palPxy ( int np, double xye[][2], double xym[][2], double coeffs[6],
              double xyp[][2], double *xrms, double *yrms, double *rrms ) {

  int i;
  double sdx2, sdy2, xp, yp, dx, dy, dx2, dy2, p;

  /* Initialize summations */
  sdx2 = 0.0;
  sdy2 = 0.0;

  /* Loop by sample */
  for (i = 0; i < np; i ++) {

    /* Transform "measured" [X,Y] to "predicted" [X,Y] */
    palXy2xy(xym[i][0], xym[i][1], coeffs, &xp, &yp);
    xyp[i][0] = xp;
    xyp[i][1] = yp;

    /* Compute residuals in X and Y, and update summations */
    dx = xye[i][0] - xp;
    dy = xye[i][1] - yp;
    dx2 = dx * dx;
    dy2 = dy * dy;
    sdx2 = sdx2 + dx2;
    sdy2 = sdy2 + dy2;

    /* Next sample */
  }

  /* Compute RMS values */
  p = DMAX(1.0, (double) np);
  *xrms = sqrt(sdx2 / p);
  *yrms = sqrt(sdy2 / p);
  *rrms = sqrt(*xrms * *xrms + *yrms * *yrms);
}
