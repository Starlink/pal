/*
*+
*  Name:
*     palXy2xy

*  Purpose:
*     Transform one [x,y] into another using a linear model of the type
*     produced by the palFitxy routine.

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palXy2xy ( double x1, double y1, double coeffs[6],
*                double *x2, double *y2)

*  Arguments:
*     x1 = double (Given)
*        x-coordinate
*     y1 = double (Given)
*        y-coordinate
*     coeffs = double[6] (Given)
*        transformation coefficients (see note)
*     x2 = double (Returned)
*        x-coordinate
*     y2 = double (Returned)
*        y-coordinate

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
*     the present routine performs the transformation:
*     ---
*        x2 = A + B * x1 + C * y1
*        y2 = D + E * x1 + F * y1
*     ---

*  See also:
*     palFitxy, palPxy, palInvf and palDcmpf

*  Authors:
*     PTW: Pat Wallace (STFC)
*     GSB: Graham Bell (EAO)

*  History:
*     1994-12-05 (PTW):
*        SLALIB implementation.
*     2018-10-23 (GSB):
*        Initial version in C.

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
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

void palXy2xy ( double x1, double y1, double coeffs[6],
                double *x2, double *y2 ) {

  *x2 = coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1;
  *y2 = coeffs[3] + coeffs[4] * x1 + coeffs[5] * y1;
}
