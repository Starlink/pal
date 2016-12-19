PAL - Positional Astronomy Library
==================================

[![DOI](https://zenodo.org/badge/5487324.svg)](https://zenodo.org/badge/latestdoi/5487324)

The PAL library is a partial re-implementation of Pat Wallace's popular SLALIB
library written in C using a Gnu GPL license and layered on top of the IAU's
SOFA library (or the BSD-licensed ERFA) where appropriate.
PAL attempts to stick to the SLA C API where
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

    ./configure --prefix=/usr/local --without-starlink
    make
    make install

The tests can be run using `make check`. Use `--prefix` to specify an install location.
Given the history of the source code as a Starlink library the default will be `/star`.

`--without-starlink` forces the configure script to forget about any Starlink
configurations. This is the safe option if you run into problems when using
a simple `--prefix` for building outside of Starlink. The configure script
will assume Starlink is not being used by looking to see if
`STARCONF_DEFAULT_PREFIX` environment variable is set. You may run into problems if
`STARCONF_DEFAULT_PREFIX` is set but you use `--without-starlink`.

Requirements
------------

Requires that either the SOFA C library or the ERFA library variant
(which has a more permissive license than SOFA) be installed.  The
`configure` script will abort if neither SOFA nor ERFA can be
found. SOFA can be obtained either from <http://www.iausofa.org> or
from an unofficial github repository (with a configure script) at
<https://github.com/Starlink/sofa/downloads>.  ERFA can be downloaded
from <https://github.com/liberfa/erfa>.

Missing Functions
-----------------

Not all SLALIB functions have been added. New routines are added to PAL as demand arises.


Language Bindings
-----------------

A Perl binding of PAL is available (<https://github.com/timj/perl-Astro-PAL>) named `Astro::PAL`
and is available from CPAN at <https://metacpan.org/module/Astro::PAL>. This is a standalone
distribution that comes with its own copies of PAL and SOFA and so can be installed directly
from the `cpan` shell.

A Python binding of PAL is available (<https://github.com/Starlink/palpy>). This is a standalone
distribution that comes with its own copies of PAL and SOFA.

The Starlink AST (<http://www.starlink.ac.uk/ast>) library now uses PAL and can be built
either with a private PAL or with an external PAL.

Documentation
-------------

The description paper for PAL is: ["_PAL: A Positional Astronomy Library_"](http://adsabs.harvard.edu/abs/2013ASPC..475..307J),
Jenness, T. & Berry, D. S., in _Astronomical Data Anaysis Software and Systems XXII_,
Friedel, D. N. (ed), ASP Conf. Ser. **475**, p307.
