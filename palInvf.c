/*
*+
*  Name:
*     palInvf

*  Purpose:
*     Invert a linear model of the type produced by the palFitxy routine.

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palInvf ( double fwds[6], double bkwds[6], int *j )

*  Arguments:
*     fwds = double[6] (Given)
*        model coefficients
*     bkwds = double[6] (Returned)
*        inverse model
*     j = int (Returned)
*        status: 0 = OK, -1 = no inverse

*  Description:
*     The models relate two sets of [x,y] coordinates as follows.
*     Naming the elements of fwds:
*     ---
*        fwds[0] = A
*        fwds[1] = B
*        fwds[2] = C
*        fwds[3] = D
*        fwds[4] = E
*        fwds[5] = F
*     ---
*     where two sets of coordinates [x1,y1] and [x2,y2] are related
*     thus:
*     ---
*        x2 = A + B * x1 + C * y1
*        y2 = D + E * x1 + F * y1
*     ---
*     the present routine generates a new set of coefficients:
*     ---
*        bkwds[0] = P
*        bkwds[1] = Q
*        bkwds[2] = R
*        bkwds[3] = S
*        bkwds[4] = T
*        bkwds[5] = U
*     ---
*     such that:
*     ---
*        x1 = P + Q * x2 + R * y2
*        y1 = S + T * x2 + U * y2
*     ---
*     Two successive calls to palInvf will thus deliver a set
*     of coefficients equal to the starting values.

*  See also:
*     palFitxy, palPxy, palXy2xy and palDcmpf

*  Authors:
*     PTW: Pat Wallace (STFC)
*     GSB: Graham Bell (EAO)

*  History:
*     1990-04-11 (PTW):
*        SLALIB implementation.
*     2004-12-26 (PTW):
*        Documentation updated.
*     2018-10-23 (GSB):
*        Initial version in C.

*  Copyright:
*     Copyright P.T.Wallace.  All rights reserved.
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

void palInvf ( double fwds[6], double bkwds[6], int *j ) {
  double a, b, c, d, e, f, det;

  a = fwds[0];
  b = fwds[1];
  c = fwds[2];
  d = fwds[3];
  e = fwds[4];
  f = fwds[5];
  det = b * f - c * e;

  if ( det != 0.0 ) {
    bkwds[0] = (c * d - a * f) / det;
    bkwds[1] = f / det;
    bkwds[2] = -c / det;
    bkwds[3] = (a * e - b * d) / det;
    bkwds[4] = -e / det;
    bkwds[5] = b / det;
    *j = 0;
  } else {
    *j = -1;
  }
}
