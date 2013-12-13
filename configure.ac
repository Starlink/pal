dnl    Process this file with autoconf to produce a configure script
AC_REVISION($Revision: 27534 $)

dnl    Initialisation: package name and version number
AC_INIT([pal],[0.3.0],[starlink@jiscmail.ac.uk])

dnl    Require autoconf-2.50 at least
AC_PREREQ([2.69])
dnl    Require Starlink automake
AM_INIT_AUTOMAKE(1.8.2-starlink)

dnl    Sanity-check: name a file in the source directory -- if this
dnl    isn't found then configure will complain
AC_CONFIG_SRCDIR([pal.h])

dnl    Include defaults for Starlink configurations
STAR_DEFAULTS

dnl    Find required versions of the programs we need for configuration
AC_PROG_CC
LT_INIT

#   If --with-pic=no is set we should honour that.
AM_CONDITIONAL(NOPIC, test x$pic_mode = xno)

dnl    copysign and isblank are a c99 feature
AC_CHECK_FUNCS(copysign)
AC_CHECK_FUNCS(isblank)

dnl    Use strlcpy if it is available
AC_CHECK_FUNCS(strlcpy)

# Need the math library
AC_CHECK_LIB([m],[sin])

dnl    We can not simply test for Starlink starutil because
dnl    when configure runs in a Starlink build starutil will not
dnl    have been built yet. If --without-starlink has been used
dnl    $STARLINK will be unset so we test that.
dnl    If we do not have Starlink we can do the test anyhow just in case

if test -n "$STARLINK"
then
  AC_SUBST( STARUTIL_LIBADD, "${libdir}/libstarutil.la" )
  AC_SUBST( SOFA_LIBADD, "${libdir}/libsofa_c.la" )
  AC_SUBST( SOFA_LDFLAGS, "" )
  AC_DEFINE( [HAVE_STAR_UTIL_H], [1], [Define to 1 if you have the <star/util.h> header file])
else
  dnl AC_CHECK_HEADERS does not search $includedir
  save_CPPFLAGS="$CPPFLAGS"
  eval CPPFLAGS=\"$CPPFLAGS -I${includdir}\"
  eval CPPFLAGS=\"$CPPFLAGS\"
  AC_CHECK_HEADERS( star/util.h )
  CPPFLAGS="$save_CPPFLAGS"

  dnl for some reason AC_CHECK_LIB does not look in the --prefix hierarchy so
  dnl $libdir is not searched.
  dnl and we use eval twice to convert $libdir -> $exec_prefix/lib -> $prefix/lib
  save_LDFLAGS="$LDFLAGS"
  eval LDFLAGS=\"$LDFLAGS -L${libdir}\"
  eval LDFLAGS=\"$LDFLAGS\"

  AC_CHECK_LIB([starutil],[star_strlcpy],
               [AC_SUBST(STARUTIL_LIBADD, "-lstarutil")],
               [AC_SUBST(STARUTIL_LIBADD, "")])

  AC_CHECK_LIB([sofa_c],[iauCal2jd],
               [AC_SUBST(SOFA_LIBADD, "-lsofa_c")],
               [AC_MSG_ERROR(Missing SOFA library. Can not continue)])
  LDFLAGS="$save_LDFLAGS"

  dnl  Ensure that we use the $prefix values
  AC_SUBST( SOFA_LDFLAGS, "-L${libdir}" )
  AC_SUBST( SOFA_CPPFLAGS, "-I${includedir}" )

  dnl  Disable document building regardless of --without-stardocs
  _star_build_docs=:
fi

dnl    Declare the build and use dependencies for this package
STAR_DECLARE_DEPENDENCIES(build, [sofa starutil])
STAR_DECLARE_DEPENDENCIES(link,  [sofa starutil])

dnl    List the sun/ssn/... numbers which document this package and
dnl    which are present as .tex files in this directory.
STAR_LATEX_DOCUMENTATION(sun267)

dnl    If you wish to configure extra files, you can add them to this
dnl    declaration.
AC_CONFIG_FILES(Makefile component.xml)
AC_CONFIG_HEADERS( config.h )

dnl    This is the bit that does the actual work
AC_OUTPUT
