/*
*+
*  Name:
*     palFitxy

*  Purpose:
*     Fit a linear model to relate two sets of [x,y] coordinates.

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palFitxy ( int itype, int np, double xye[][2], double xym[][2],
*                double coeffs[6], int *j )

*  Arguments:
*     itype = int (Given)
*        type of model: 4 or 6 (note 1)
*     np = int (Given)
*        number of samples (note 2)
*     xye = double[np][2] (Given)
*        expected [x,y] for each sample
*     xym = double[np][2] (Given)
*        measured [x,y] for each sample
*     coeffs = double[6] (Returned)
*        coefficients of model (note 3)
*     j = int * (Returned)
*        status:
*        -  0 = OK
*        - -1 = illegal itype
*        - -2 = insufficient data
*        - -3 = no solution

*  Description:
*     Fits a linear model to relate two sets of [X,Y] coordinates.

*  Notes:
*     1)  itype, which must be either 4 or 6, selects the type of model
*         fitted.  Both allowed itype values produce a model coeffs which
*         consists of six coefficients, namely the zero points and, for
*         each of xe and ye, the coefficient of xm and ym.  For itype=6,
*         all six coefficients are independent, modelling squash and shear
*         as well as origin, scale, and orientation.  However, itype=4
*         selects the "solid body rotation" option;  the model coeffs
*         still consists of the same six coefficients, but now two of
*         them are used twice (appropriately signed).  Origin, scale
*         and orientation are still modelled, but not squash or shear -
*         the units of x and y have to be the same.
*
*     2)  For itype=4, np must be at least 2.  For itype=6, np must be at
*         least 3.
*
*     3)  The model is returned in the array coeffs.  Naming the
*         elements of coeffs as follows:
*         ---
*            coeffs[0] = A
*            coeffs[1] = B
*            coeffs[2] = C
*            coeffs[3] = D
*            coeffs[4] = E
*            coeffs[5] = F
*         ---
*         the model is:
*         ---
*            xe = A + B * xm + C * ym
*            ye = D + E * xm + F * ym
*         ---
*         For the "solid body rotation" option (itype=4), the
*         magnitudes of B and F, and of C and E, are equal.  The
*         signs of these coefficients depend on whether there is a
*         sign reversal between xe,ye and xm,ym;  fits are performed
*         with and without a sign reversal and the best one chosen.
*
*     4)  Error status values j=-1 and -2 leave coeffs unchanged;
*         if j=-3 coeffs may have been changed.

*  See also:
*      palPxy, palInvf, palXy2xy and palDcmpf

*  Authors:
*     PTW: Pat Wallace (STFC)
*     GSB: Graham Bell (EAO)

*  History:
*     2001-11-30 (PTW):
*        SLALIB implementation.
*     2005-09-08 (PTW):
*        Fix compiler uninitialised warnings.
*     2018-10-23 (GSB):
*        Initial version in C.

*  Copyright:
*     Copyright (C) 2005 P.T.Wallace. All rights reserved.
*     Copyright (C) 2018 East Asian Observatory.
*
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

#include "pal.h"

void palFitxy ( int itype, int np, double xye[][2], double xym[][2],
                double coeffs[6], int *j) {

  int i, jstat, iw[4], nsol;
  double a, b, c, d, aold, bold, cold, dold, sold,
    p, sxe, sxexm, sxeym, sye, syeym, syexm, sxm,
    sym, sxmxm, sxmym, symym, xe, ye,
    xm, ym, v[4], dm3[3][3], dm4[4][4], det,
    sgn, sxxyy, sxyyx, sx2y2, sdr2, xr, yr;

  /* Preset the status */
  *j = 0;

  /* Variable initializations to avoid compiler warnings */
  a = 0.0;
  b = 0.0;
  c = 0.0;
  d = 0.0;
  aold = 0.0;
  bold = 0.0;
  cold = 0.0;
  dold = 0.0;
  sold = 0.0;

  /* Float the number of samples */
  p = (double) np;

  /* Check ITYPE */
  if (itype == 6) {
    /* Six-coefficient linear model */

    /* Check enough samples */
    if (np >= 3) {

      /* Form summations */
      sxe = 0.0;
      sxexm = 0.0;
      sxeym = 0.0;
      sye = 0.0;
      syeym = 0.0;
      syexm = 0.0;
      sxm = 0.0;
      sym = 0.0;
      sxmxm = 0.0;
      sxmym = 0.0;
      symym = 0.0;

      for (i = 0; i < np; i++) {
        xe = xye[i][0];
        ye = xye[i][1];
        xm = xym[i][0];
        ym = xym[i][1];
        sxe = sxe + xe;
        sxexm = sxexm + xe * xm;
        sxeym = sxeym + xe * ym;
        sye = sye + ye;
        syeym = syeym + ye * ym;
        syexm = syexm + ye * xm;
        sxm = sxm + xm;
        sym = sym + ym;
        sxmxm = sxmxm + xm * xm;
        sxmym = sxmym + xm * ym;
        symym = symym + ym * ym;
      }

      /* Solve for A, B, C in  xe = A + B * xm + C * ym */
      v[0] = sxe;
      v[1] = sxexm;
      v[2] = sxeym;
      dm3[0][0] = p;
      dm3[0][1] = sxm;
      dm3[0][2] = sym;
      dm3[1][0] = sxm;
      dm3[1][1] = sxmxm;
      dm3[1][2] = sxmym;
      dm3[2][0] = sym;
      dm3[2][1] = sxmym;
      dm3[2][2] = symym;

      palDmat(3, *dm3, v, &det, &jstat, iw);

      if (jstat == 0) {
         for (i = 0; i < 3; i ++) {
            coeffs[i] = v[i];
         }

         /* Solve for D, E, F in  ye = D + E * xm + F * ym */
         v[0] = sye;
         v[1] = syexm;
         v[2] = syeym;

         palDmxv(dm3, v, coeffs + 3);
      } else {
        /* No 6-coefficient solution possible */
        *j = -3;
      }
    } else {
      /* Insufficient data for 6-coefficient fit */
      *j = -2;
    }
  } else if (itype == 4) {
    /* Four-coefficient solid body rotation model */

    /* Check enough samples */
    if (np >= 2) {

      /* Try two solutions, first without then with flip in X */
      for (nsol = 0; nsol < 2; nsol ++) {
        if (nsol == 0) {
          sgn = 1.0;
        } else {
          sgn = -1.0;
        }

        /* Form summations */
        sxe = 0.0;
        sxxyy = 0.0;
        sxyyx = 0.0;
        sye = 0.0;
        sxm = 0.0;
        sym = 0.0;
        sx2y2 = 0.0;

        for (i = 0; i < np; i ++) {
          xe = xye[i][0] * sgn;
          ye = xye[i][1];
          xm = xym[i][0];
          ym = xym[i][1];
          sxe = sxe + xe;
          sxxyy = sxxyy + xe * xm + ye * ym;
          sxyyx = sxyyx + xe * ym - ye * xm;
          sye = sye + ye;
          sxm = sxm + xm;
          sym = sym + ym;
          sx2y2 = sx2y2 + xm * xm + ym * ym;
        }

        /* Solve for A, B, C, D in:  +/- xe = A + B * xm - C * ym
                                       + ye = D + C * xm + B * ym */
        v[0] = sxe;
        v[1] = sxxyy;
        v[2] = sxyyx;
        v[3] = sye;
        dm4[0][0] = p;
        dm4[0][1] = sxm;
        dm4[0][2] = -sym;
        dm4[0][3] = 0.0;
        dm4[1][0] = sxm;
        dm4[1][1] = sx2y2;
        dm4[1][2] = 0.0;
        dm4[1][3] = sym;
        dm4[2][0] = sym;
        dm4[2][1] = 0.0;
        dm4[2][2] = -sx2y2;
        dm4[2][3] = -sxm;
        dm4[3][0] = 0.0;
        dm4[3][1] = sym;
        dm4[3][2] = sxm;
        dm4[3][3] = p;

        palDmat(4, *dm4, v, &det, &jstat, iw);

        if (jstat == 0) {
          a = v[0];
          b = v[1];
          c = v[2];
          d = v[3];

          /* Determine sum of radial errors squared */
          sdr2 = 0.0;
          for (i = 0; i < np; i ++) {
            xm = xym[i][0];
            ym = xym[i][1];
            xr = a + b * xm - c * ym - xye[i][0] * sgn;
            yr = d + c * xm + b * ym - xye[i][1];
            sdr2 = sdr2 + xr * xr + yr * yr;
          }

        } else {
          /* Singular: set flag */
          sdr2 = -1.0;
        }

        /* If first pass and non-singular, save variables */
        if (nsol == 0 && jstat == 0) {
          aold = a;
          bold = b;
          cold = c;
          dold = d;
          sold = sdr2;
        }
      }

      /* Pick the best of the two solutions */
      if (sold >= 0.0 && (sold <= sdr2 || np == 2)) {
        coeffs[0] = aold;
        coeffs[1] = bold;
        coeffs[2] = -cold;
        coeffs[3] = dold;
        coeffs[4] = cold;
        coeffs[5] = bold;
      } else if (jstat == 0) {
        coeffs[0] = -a;
        coeffs[1] = -b;
        coeffs[2] = c;
        coeffs[3] = d;
        coeffs[4] = c;
        coeffs[5] = b;
      } else {
        /* No 4-coefficient fit possible */
        *j = -3;
      }
    } else {
      /* Insufficient data for 4-coefficient fit */
      *j = -2;
    }
  } else {
    /* Illegal itype - not 4 or 6 */
    *j = -1;
  }
}
