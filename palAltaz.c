/*
*+
*  Name:
*     palAltaz

*  Purpose:
*     Positions, velocities and accelerations for an altazimuth telescope mount

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palAltaz ( double ha, double dec, double phi,
*                double *az, double *azd, double *azdd,
*                double *el, double *eld, double *eldd,
*                double *pa, double *pad, double *padd );

*  Arguments:
*     ha = double (Given)
*        Hour angle (radians)
*     dec = double (Given)
*        Declination (radians)
*     phi = double (Given)
*        Observatory latitude (radians)
*     az = double * (Returned)
         Azimuth (radians)
*     azd = double * (Returned)
*        Azimuth velocity (radians per radian of HA)
*     azdd = double * (Returned)
*        Azimuth acceleration (radians per radian of HA squared)
*     el = double * (Returned)
*        Elevation (radians)
*     eld = double * (Returned)
*        Elevation velocity (radians per radian of HA)
*     eldd = double * (Returned)
*        Elevation acceleration (radians per radian of HA squared)
*     pa = double * (Returned)
*        Parallactic angle (radians)
*     pad = double * (Returned)
*        Parallactic angle velocity (radians per radian of HA)
*     padd = double * (Returned)
*        Parallactic angle acceleration (radians per radian of HA squared)


*  Description:
*     Positions, velocities and accelerations for an altazimuth
*     telescope mount.

*  Authors:
*     PTW: P. T. Wallace
*     TIMJ: Tim Jenness (Cornell)
*     {enter_new_authors_here}

*  Notes:
*     - Natural units are used throughout.  HA, DEC, PHI, AZ, EL
*       and ZD are in radians.  The velocities and accelerations
*       assume constant declination and constant rate of change of
*       hour angle (as for tracking a star);  the units of AZD, ELD
*       and PAD are radians per radian of HA, while the units of AZDD,
*       ELDD and PADD are radians per radian of HA squared.  To
*       convert into practical degree- and second-based units:
*
*         angles * 360/2pi -> degrees
*         velocities * (2pi/86400)*(360/2pi) -> degree/sec
*         accelerations * ((2pi/86400)**2)*(360/2pi) -> degree/sec/sec
*
*       Note that the seconds here are sidereal rather than SI.  One
*       sidereal second is about 0.99727 SI seconds.
*
*       The velocity and acceleration factors assume the sidereal
*       tracking case.  Their respective numerical values are (exactly)
*       1/240 and (approximately) 1/3300236.9.
*
*     - Azimuth is returned in the range 0-2pi;  north is zero,
*       and east is +pi/2.  Elevation and parallactic angle are
*       returned in the range +/-pi.  Parallactic angle is +ve for
*       a star west of the meridian and is the angle NP-star-zenith.
*
*     - The latitude is geodetic as opposed to geocentric.  The
*       hour angle and declination are topocentric.  Refraction and
*       deficiencies in the telescope mounting are ignored.  The
*       purpose of the routine is to give the general form of the
*       quantities.  The details of a real telescope could profoundly
*       change the results, especially close to the zenith.
*
*     - No range checking of arguments is carried out.
*
*     - In applications which involve many such calculations, rather
*       than calling the present routine it will be more efficient to
*       use inline code, having previously computed fixed terms such
*       as sine and cosine of latitude, and (for tracking a star)
*       sine and cosine of declination.

*  History:
*     2014-09-30 (TIMJ):
*        Initial version. Ported from Fortran SLA
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2004 P.T. Wallace
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
#include "palmac.h"

void
palAltaz ( double ha, double dec, double phi,
           double *az, double *azd, double *azdd,
           double *el, double *eld, double *eldd,
           double *pa, double *pad, double *padd ) {

  const double TINY = 1E-30;

  double sh,ch,sd,cd,sp,cp,chcd,sdcp,x,y,z,rsq,r,a,e,c,s,
    q,qd,ad,ed,edr,add,edd,qdd;


  /*  Useful functions */
  sh=sin(ha);
  ch=cos(ha);
  sd=sin(dec);
  cd=cos(dec);
  sp=sin(phi);
  cp=cos(phi);
  chcd=ch*cd;
  sdcp=sd*cp;
  x=-chcd*sp+sdcp;
  y=-sh*cd;
  z=chcd*cp+sd*sp;
  rsq=x*x+y*y;
  r=sqrt(rsq);

  /*  Azimuth and elevation */
  if (rsq == 0.0) {
    a=0.0;
  } else {
    a=atan2(y,x);
  }
  if (a < 0.0) a += PAL__D2PI;
  e=atan2(z,r);

  /*  Parallactic angle */
  c=cd*sp-ch*sdcp;
  s=sh*cp;
  if (c*c+s*s > 0) {
    q=atan2(s,c);
  } else {
    q= PAL__DPI - ha;
  }

  /*  Velocities and accelerations (clamped at zenith/nadir) */
  if (rsq < TINY) {
    rsq=TINY;
    r=sqrt(rsq);
  }
  qd=-x*cp/rsq;
  ad=sp+z*qd;
  ed=cp*y/r;
  edr=ed/r;
  add=edr*(z*sp+(2.0-rsq)*qd);
  edd=-r*qd*ad;
  qdd=edr*(sp+2.0*z*qd);

  /*  Results */
  *az=a;
  *azd=ad;
  *azdd=add;
  *el=e;
  *eld=ed;
  *eldd=edd;
  *pa=q;
  *pad=qd;
  *padd=qdd;

}
