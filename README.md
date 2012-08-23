PAL - Positional Astronomy Library
==================================

The PAL library is a partial re-implementation of Pat Wallace's popular SLALIB
library written in C using a Gnu GPL license and layered on top of the IAU's
SOFA library where appropriate. PAL attempts to stick to the SLA C API where
possible although `palObs()` has a more C-like API than the equivalent
`slaObs()` function. In most cases it is enough to simply change the function
prefix of a routine in order to link against PAL rather than SLALIB. Routines
calling SOFA use modern nutation and precession models so will return slightly
different answers than native SLALIB. PAL functions not available in SOFA were
ported from the Fortran version of SLALIB that ships as part of the Starlink
software and uses a GPL licence.

See `pal.news` for release notes.

Building
--------

A simple `configure` script is provided:

    ./configure --prefix=/usr/local
    make
    make install

The tests can be run using `make check`. Use `--prefix` to specify an install location.
Given the history of the source code as a Starlink libary the default will be `/star`.

Requirements
------------

Requires that the SOFA C library is installed. The `configure` script will abort if SOFA
can not be found. SOFA can be obtained either from <http://www.iausofa.org> or from
an unofficial github repository (with a configure script) at <https://github.com/Starlink/sofa/downloads>.

Missing Functions
-----------------

Not all SLALIB functions have been added. New routines are added to PAL as demand arises.


Language Bindings
-----------------

A Perl binding of PAL is available (<https://github.com/timj/perl-Astro-PAL>) named `Astro::PAL`
and is available from CPAN at <https://metacpan.org/module/Astro::PAL>. This is a standalone
distribution that comes with its own copies of PAL and SOFA and so can be installed directly
from the `cpan` shell.

The Starlink AST (<http://www.starlink.ac.uk/ast>) library now uses PAL and can be built
either with a private PAL or with an external PAL.