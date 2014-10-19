/*
*+
*  Name:
*     palPolmo

*  Purpose:
*     Correct for polar motion

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palPolmo ( double elongm, double phim, double xp, double yp,
*                double *elong, double *phi, double *daz );

*  Arguments:
*     elongm = double (Given)
*        Mean logitude of the observer (radians, east +ve)
*     phim = double (Given)
*        Mean geodetic latitude of the observer (radians)
*     xp = double (Given)
*        Polar motion x-coordinate (radians)
*     yp = double (Given)
*        Polar motion y-coordinate (radians)
*     elong = double * (Returned)
*        True longitude of the observer (radians, east +ve)
*     phi = double * (Returned)
*        True geodetic latitude of the observer (radians)
*     daz = double * (Returned)
*        Azimuth correction (terrestrial-celestial, radians)

*  Description:
*     Polar motion:  correct site longitude and latitude for polar
*     motion and calculate azimuth difference between celestial and
*     terrestrial poles.

*  Authors:
*     PTW: Patrick Wallace (STFC)
*     TIMJ: Tim Jenness (Cornell)
*     {enter_new_authors_here}

*  Notes:
*     - "Mean" longitude and latitude are the (fixed) values for the
*       site's location with respect to the IERS terrestrial reference
*       frame;  the latitude is geodetic.  TAKE CARE WITH THE LONGITUDE
*       SIGN CONVENTION.  The longitudes used by the present routine
*       are east-positive, in accordance with geographical convention
*       (and right-handed).  In particular, note that the longitudes
*       returned by the sla_OBS routine are west-positive, following
*       astronomical usage, and must be reversed in sign before use in
*       the present routine.
*
*     - XP and YP are the (changing) coordinates of the Celestial
*       Ephemeris Pole with respect to the IERS Reference Pole.
*       XP is positive along the meridian at longitude 0 degrees,
*       and YP is positive along the meridian at longitude
*       270 degrees (i.e. 90 degrees west).  Values for XP,YP can
*       be obtained from IERS circulars and equivalent publications;
*       the maximum amplitude observed so far is about 0.3 arcseconds.
*
*     - "True" longitude and latitude are the (moving) values for
*       the site's location with respect to the celestial ephemeris
*       pole and the meridian which corresponds to the Greenwich
*       apparent sidereal time.  The true longitude and latitude
*       link the terrestrial coordinates with the standard celestial
*       models (for precession, nutation, sidereal time etc).
*
*     - The azimuths produced by sla_AOP and sla_AOPQK are with
*       respect to due north as defined by the Celestial Ephemeris
*       Pole, and can therefore be called "celestial azimuths".
*       However, a telescope fixed to the Earth measures azimuth
*       essentially with respect to due north as defined by the
*       IERS Reference Pole, and can therefore be called "terrestrial
*       azimuth".  Uncorrected, this would manifest itself as a
*       changing "azimuth zero-point error".  The value DAZ is the
*       correction to be added to a celestial azimuth to produce
*       a terrestrial azimuth.
*
*     - The present routine is rigorous.  For most practical
*       purposes, the following simplified formulae provide an
*       adequate approximation:
*
*       elong = elongm+xp*cos(elongm)-yp*sin(elongm)
*       phi   = phim+(xp*sin(elongm)+yp*cos(elongm))*tan(phim)
*       daz   = -sqrt(xp*xp+yp*yp)*cos(elongm-atan2(xp,yp))/cos(phim)
*
*       An alternative formulation for DAZ is:
*
*       x = cos(elongm)*cos(phim)
*       y = sin(elongm)*cos(phim)
*       daz = atan2(-x*yp-y*xp,x*x+y*y)
*
*     - Reference:  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement
*                   to the Astronomical Almanac", ISBN 0-935702-68-7,
*                   sections 3.27, 4.25, 4.52.

*  History:
*     2000-11-30 (PTW):
*        SLALIB implementation.
*     2014-10-18 (TIMJ):
*        Initial version in C.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2000 Rutherford Appleton Laboratory.
*     Copyright (C) 2014 Cornell University
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

#include <math.h>

#include "pal.h"

void palPolmo ( double elongm, double phim, double xp, double yp,
                double *elong, double *phi, double *daz ) {

  double  sel,cel,sph,cph,xm,ym,zm,xnm,ynm,znm,
    sxp,cxp,syp,cyp,zw,xt,yt,zt,xnt,ynt;

  /*  Site mean longitude and mean geodetic latitude as a Cartesian vector */
  sel=sin(elongm);
  cel=cos(elongm);
  sph=sin(phim);
  cph=cos(phim);

  xm=cel*cph;
  ym=sel*cph;
  zm=sph;

  /*  Rotate site vector by polar motion, Y-component then X-component */
  sxp=sin(xp);
  cxp=cos(xp);
  syp=sin(yp);
  cyp=cos(yp);

  zw=(-ym*syp+zm*cyp);

  xt=xm*cxp-zw*sxp;
  yt=ym*cyp+zm*syp;
  zt=xm*sxp+zw*cxp;

  /*  Rotate also the geocentric direction of the terrestrial pole (0,0,1) */
  xnm=-sxp*cyp;
  ynm=syp;
  znm=cxp*cyp;

  cph=sqrt(xt*xt+yt*yt);
  if (cph == 0.0) xt=1.0;
  sel=yt/cph;
  cel=xt/cph;

  /*  Return true longitude and true geodetic latitude of site */
  if (xt != 0.0 || yt != 0.0) {
    *elong=atan2(yt,xt);
  } else {
    *elong=0.0;
  }
  *phi=atan2(zt,cph);

  /*  Return current azimuth of terrestrial pole seen from site position */
  xnt=(xnm*cel+ynm*sel)*zt-znm*cph;
  ynt=-xnm*sel+ynm*cel;
  if (xnt != 0.0 || ynt != 0.0) {
    *daz=atan2(-ynt,-xnt);
  } else {
    *daz=0.0;
  }

}
