/*
*+
*  Name:
*     palVers

*  Purpose:
*     Obtain PAL version number

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     int palVers( char *verstring, size_t verlen );

*  Arguments:
*     verstring = char * (Returned)
*        Buffer to receive version string of the form "A.B.C". Can be NULL.
*     verlen = size_t (Given)
*        Allocated size of "verstring" including nul. Version string
*        will be truncated if it does not fit in buffer.

*  Returned Value:
*     vernum = int (Returned)
*        Version number as an integer.

*  Description:
*     Retrieve the PAL version number as a string in the form "A.B.C"
*     and as an integer (major*1e6+minor*1e3+release).

*  Authors:
*     TIMJ: Tim Jenness (Cornell)
*     {enter_new_authors_here}

*  Notes:
*     - Note that this API does not match the slaVers API.

*  History:
*     2014-08-27 (TIMJ):
*        Initial version
*     {enter_further_changes_here}

*  Copyright:
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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <string.h>

/* This version is just a straight copy without putting ellipsis on the end. */
static void my__strlcpy( char * dest, const char * src, size_t size ) {
# if HAVE_STRLCPY
  strlcpy( dest, src, size );
# else
  strncpy( dest, src, size );
  dest[size-1] = '\0';
# endif
}

int
palVers( char *verstring, size_t verlen ) {
  if (verstring) my__strlcpy( verstring, PACKAGE_VERSION, verlen );
  return PACKAGE_VERSION_INTEGER;
}
