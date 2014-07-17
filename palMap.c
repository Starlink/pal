/*
*+
*  Name:
*     palMap

*  Purpose:
*     Convert star RA,Dec from mean place to geocentric apparent

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     void palMap( double rm, double dm, double pr, double pd,
*                  double px, double rv, double eq, double date,
*                  double *ra, double *da );

*  Arguments:
*     rm = double (Given)
*        Mean RA (radians)
*     dm = double (Given)
*        Mean declination (radians)
*     pr = double (Given)
*        RA proper motion, changes per Julian year (radians)
*     pd = double (Given)
*        Dec proper motion, changes per Julian year (radians)
*     px = double (Given)
*        Parallax (arcsec)
*     rv = double (Given)
*        Radial velocity (km/s, +ve if receding)
*     eq = double (Given)
*        Epoch and equinox of star data (Julian)
*     date = double (Given)
*        TDB for apparent place (JD-2400000.5)
*     ra = double * (Returned)
*        Apparent RA (radians)
*     dec = double * (Returned)
*        Apparent dec (radians)

*  Description:
*     Convert star RA,Dec from mean place to geocentric apparent.

*  Authors:
*     PTW: Patrick T. Wallace
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - Calls palMappa and palMapqk

*     - The reference systems and timescales used are IAU 2006.

*     - EQ is the Julian epoch specifying both the reference frame and
*       the epoch of the position - usually 2000.  For positions where
*       the epoch and equinox are different, use the routine palPm to
*       apply proper motion corrections before using this routine.
*
*     - The distinction between the required TDB and TT is always
*       negligible.  Moreover, for all but the most critical
*       applications UTC is adequate.
*
*     - The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt.
*
*     - This routine may be wasteful for some applications because it
*       recomputes the Earth position/velocity and the precession-
*       nutation matrix each time, and because it allows for parallax
*       and proper motion.  Where multiple transformations are to be
*       carried out for one epoch, a faster method is to call the
*       palMappa routine once and then either the palMapqk routine
*       (which includes parallax and proper motion) or palMapqkz (which
*       assumes zero parallax and proper motion).
*
*     - The accuracy is sub-milliarcsecond, limited by the
*       precession-nutation model (see palPrenut for details).
*
*     - The accuracy is further limited by the routine palEvp, called
*       by palMappa, which computes the Earth position and velocity.
*       See iauEpv00 for details on that calculation.

*  History:
*     2012-03-01 (TIMJ):
*        Initial version
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 2001 Rutherford Appleton Laboratory
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
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
*    USA.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#include "pal.h"

void palMap( double rm, double dm, double pr, double pd,
             double px, double rv, double eq, double date,
             double *ra, double *da ) {

  double amprms[21];

  /* Star independent parameters */
  palMappa( eq, date, amprms );

  /* Mean to apparent */
  palMapqk( rm, dm, pr, pd, px, rv, amprms, ra, da );

}
