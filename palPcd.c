/*
*+
*  Name:
*     palPcd

*  Purpose:
*     Apply pincushion/barrel distortion to a tangent-plane [x,y]

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palPcd( double disco, double * x, double * y );

*  Arguments:
*     disco = double (Given)
*        Pincushion/barrel distortion coefficient.
*     x = double * (Given & Returned)
*        On input the tangent-plane X coordinate, on output
*        the distorted X coordinate.
*     y = double * (Given & Returned)
*        On input the tangent-plane Y coordinate, on output
*        the distorted Y coordinate.

*  Description:
*     Applies pincushion and barrel distortion to a tangent
*     plane coordinate.

*  Authors:
*     PTW: Pat Wallace (RAL)
*     TIMJ: Tim Jenness (Cornell)
*     {enter_new_authors_here}

*  Notes:
*     - The distortion is of the form RP = R*(1 + C*R**2), where R is
*       the radial distance from the tangent point, C is the DISCO
*       argument, and RP is the radial distance in the presence of
*       the distortion.
*
*     - For pincushion distortion, C is +ve;  for barrel distortion,
*       C is -ve.
*
*     - For X,Y in units of one projection radius (in the case of
*       a photographic plate, the focal length), the following
*       DISCO values apply:
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
*  See Also:
*     - There is a companion routine, palUnpcd, which performs the
*       inverse operation.

*  History:
*     2000-09-03 (PTW):
*        SLALIB implementation.
*     2015-01-01 (TIMJ):
*        Initial version. Ported from Fortran.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2000 Rutherford Appleton Laboratory.
*     Copyright (C) 2015 Cornell University
*     All Rights Reserved.

*  Licence:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 3 of the License, or
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

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#include "pal.h"

void palPcd( double disco, double *x, double *y ) {
  double f;

  f = 1.0 + disco * ( (*x) * (*x) + (*y) * (*y) );
  *x *= f;
  *y *= f;
}

