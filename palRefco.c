/*
*+
*  Name:
*     palRefco

*  Purpose:
*     Determine constants in atmospheric refraction model

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palRefco ( double hm, double tdk, double pmb, double rh,
*                     double wl, double phi, double tlr, double eps,
*                     double *refa, double *refb );

*  Arguments:
*     hm = double (Given)
*        Height of the observer above sea level (metre)
*     tdk = double (Given)
*        Ambient temperature at the observer (K)
*     pmb = double (Given)
*        Pressure at the observer (millibar)
*     rh = double (Given)
*        Relative humidity at the observer (range 0-1)
*     wl = double (Given)
*        Effective wavelength of the source (micrometre)
*     phi = double (Given)
*        Latitude of the observer (radian, astronomical)
*     tlr = double (Given)
*        Temperature lapse rate in the troposphere (K/metre)
*     eps = double (Given)
*        Precision required to terminate iteration (radian)
*     refa = double * (Returned)
*        tan Z coefficient (radian)
*     refb = double * (Returned)
*        tan**3 Z coefficient (radian)

*  Description:
*     Determine the constants A and B in the atmospheric refraction
*     model dZ = A tan Z + B tan**3 Z.
*
*     Z is the "observed" zenith distance (i.e. affected by refraction)
*     and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
*     zenith distance.

*  Authors:
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - Typical values for the TLR and EPS arguments might be 0.0065D0 and
*     1D-10 respectively.
*
*     - The radio refraction is chosen by specifying WL > 100 micrometres.
*
*     - The routine is a slower but more accurate alternative to the
*     palRefcoq routine.  The constants it produces give perfect
*     agreement with palRefro at zenith distances arctan(1) (45 deg)
*     and arctan(4) (about 76 deg).  It achieves 0.5 arcsec accuracy
*     for ZD < 80 deg, 0.01 arcsec accuracy for ZD < 60 deg, and
*     0.001 arcsec accuracy for ZD < 45 deg.

*  History:
*     2012-08-24 (TIMJ):
*        Initial version. A direct copy of the Fortran SLA implementation.
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
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

void palRefco ( double hm, double tdk, double pmb, double rh,
                double wl, double phi, double tlr, double eps,
                double *refa, double *refb ) {

  double r1, r2;

  /*  Sample zenith distances: arctan(1) and arctan(4) */
  const double ATN1 = 0.7853981633974483;
  const double ATN4 = 1.325817663668033;

  /*  Determine refraction for the two sample zenith distances */
  palRefro(ATN1,hm,tdk,pmb,rh,wl,phi,tlr,eps,&r1);
  palRefro(ATN4,hm,tdk,pmb,rh,wl,phi,tlr,eps,&r2);

  /*  Solve for refraction constants */
  *refa = (64.0*r1-r2)/60.0;
  *refb = (r2-4.0*r1)/60.0;

}
