dnl    Process this file with autoconf to produce a configure script
AC_REVISION($Revision: 27534 $)

dnl    Initialisation: package name and version number
AC_INIT([pal],[0.9.3],[starlink@jiscmail.ac.uk])
AC_CONFIG_AUX_DIR([build-aux])

dnl    Require autoconf-2.50 at least
AC_PREREQ([2.69])
dnl    Require Starlink automake
AM_INIT_AUTOMAKE(1.8.2-starlink)

dnl    Sanity-check: name a file in the source directory -- if this
dnl    isn't found then configure will complain
AC_CONFIG_SRCDIR([pal.h])

dnl    Include defaults for Starlink configurations
STAR_DEFAULTS

dnl    Would like the version number as an integer
AC_DEFINE_UNQUOTED([PACKAGE_VERSION_INTEGER], $PACKAGE_VERSION_INTEGER,
                    [Integer version number, in the form major*1e6+minor*1e3+release])

dnl    Find required versions of the programs we need for configuration
AC_PROG_CC
LT_INIT

#   If --with-pic=no is set we should honour that.
AM_CONDITIONAL(NOPIC, test x$pic_mode = xno)

dnl    copysign and isblank are a c99 feature
AC_CHECK_FUNCS(copysign)
AC_CHECK_FUNCS(isblank)

dnl    Use strlcpy if it is available
AC_SEARCH_LIBS([strlcpy], [bsd])
AS_IF([test "x$ac_cv_search_strlcpy" = "x-lbsd"], [AC_CHECK_HEADERS([bsd/string.h])])
AC_CHECK_FUNCS([strlcpy])

# Need the math library
AC_CHECK_LIB([m],[sin])

dnl    We can not simply test for Starlink starutil because
dnl    when configure runs in a Starlink build starutil will not
dnl    have been built yet. If --without-starlink has been used
dnl    $STARLINK will be unset but to play safe we also check STARCONF_DEFAULT_PREFIX
dnl    If we do not have Starlink we can do the test anyhow just in case

if test -n "$STARCONF_DEFAULT_PREFIX" -a -n "$STARLINK"
then
  AC_MSG_NOTICE([Assuming a Starlink environment])
  AC_SUBST( STARUTIL_LIBADD, "${libdir}/libstarutil.la" )
  AC_SUBST( ERFA_LIBADD, "${libdir}/liberfa.la" )
  AC_SUBST( ERFA_LDFLAGS, "" )
  AC_DEFINE( [HAVE_STAR_UTIL_H], [1], [Define to 1 if you have the <star/util.h> header file])
else

  #   Allow ERFA/SOFA location to be specified using --with-erfa=$ERFA_DIR
  #   Assumes that the value supplied here is the root and lib and include directories
  #   are below it. --with-erfa=no or --without-erfa will result in ERFA being
  #   located in $PREFIX tree. This option is only effective if Starlink is not
  #   active.
  AC_ARG_WITH(erfa,
              AS_HELP_STRING([--with-erfa],
                             [Location of ERFA/SOFA tree]),
              [if test -z "$withval" -o "$withval" = yes; then
                   unset ERFA_DIR
               elif test "X$withval" = Xno; then
                   unset ERFA_DIR
               elif test -d "$withval"; then
                   ERFA_DIR="$withval"
               else
                   unset ERFA_DIR
                   AC_MSG_WARN([--with-erfa given nonexistent directory; ignored: using default instead])
               fi])
  if test -n "$ERFA_DIR"; then
      AC_MSG_NOTICE([ERFA/SOFA tree located at $ERFA_DIR])
      erfa_includedir="${ERFA_DIR}/include"
      erfa_libdir="${ERFA_DIR}/lib"
  else
      ERFA_DIR=${prefix}
      erfa_includedir=${includedir}
      erfa_libdir=${libdir}
      AC_MSG_NOTICE([Looking for ERFA/SOFA in default location of $ERFA_DIR])
  fi

  dnl AC_CHECK_HEADERS does not search $includedir
  save_CPPFLAGS="$CPPFLAGS"
  eval CPPFLAGS=\"$CPPFLAGS -I${includedir} -I${erfa_includedir}\"
  eval CPPFLAGS=\"$CPPFLAGS\"
  AC_CHECK_HEADERS( star/util.h )
  CPPFLAGS="$save_CPPFLAGS"

  dnl for some reason AC_CHECK_LIB does not look in the --prefix hierarchy so
  dnl $libdir is not searched.
  dnl and we use eval twice to convert $libdir -> $exec_prefix/lib -> $prefix/lib
  save_LDFLAGS="$LDFLAGS"
  eval LDFLAGS=\"$LDFLAGS -L${libdir} -L${erfa_libdir}\"
  eval LDFLAGS=\"$LDFLAGS\"

  AC_CHECK_LIB([starutil],[star_strlcpy],
               [AC_SUBST(STARUTIL_LIBADD, "-lstarutil")],
               [AC_SUBST(STARUTIL_LIBADD, "")])

  AC_CHECK_LIB([erfa],[eraCal2jd],
               [AC_SUBST(ERFA_LIBADD, "-lerfa")],
               [
                AC_CHECK_LIB([sofa_c],[iauCal2jd],
                             [AC_SUBST(ERFA_LIBADD, "-lsofa_c")
                              AC_DEFINE([HAVE_SOFA_H],[1],"Build with SOFA library")],
                             [AC_MSG_ERROR(Neither ERFA nor SOFA library located. Can not continue)])
               ])
  LDFLAGS="$save_LDFLAGS"

  dnl  Ensure that we use the $prefix values and the ERFA values
  AC_SUBST( ERFA_LDFLAGS, "-L${libdir} -L${erfa_libdir}" )
  AC_SUBST( ERFA_CPPFLAGS, "-I${includedir} -I${erfa_includedir}" )

  dnl  Disable document building regardless of --without-stardocs
  _star_build_docs=:
fi

dnl    Declare the build and use dependencies for this package
STAR_DECLARE_DEPENDENCIES(build, [erfa starutil])
STAR_DECLARE_DEPENDENCIES(link,  [erfa starutil])

dnl    List the sun/ssn/... numbers which document this package and
dnl    which are present as .tex files in this directory.
STAR_LATEX_DOCUMENTATION(sun267)

dnl    If you wish to configure extra files, you can add them to this
dnl    declaration.
AC_CONFIG_FILES(Makefile component.xml)
AC_CONFIG_HEADERS( config.h )

dnl    This is the bit that does the actual work
AC_OUTPUT
