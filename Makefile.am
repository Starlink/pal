## Process this file with automake to produce Makefile.in

lib_LTLIBRARIES = libpal.la


# Include palOne2One.c separately since it is a combo file
libpal_la_SOURCES = $(PUBLIC_C_FILES) palOne2One.c $(PRIVATE_C_FILES)

# If we are using a non-standard location
libpal_la_CPPFLAGS = $(ERFA_CPPFLAGS)
libpal_la_LDFLAGS = $(ERFA_LDFLAGS)

# Force a link against ERFA and, optionally, starutil
libpal_la_LIBADD = $(ERFA_LIBADD) $(STARUTIL_LIBADD)

# Misc files
dist_starnews_DATA = pal.news
dist_pkgdata_DATA = COPYING COPYING.LESSER README.md

# Make all library code position independent. This is handy for creating
# shareable libraries from the static ones (Java JNI libraries).
if !NOPIC
libpal_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
endif

# install pal as "star/pal.h"
cincludedir = $(includedir)/star
cinclude_HEADERS = pal.h palmac.h

noinst_HEADERS = $(PRIVATE_INCLUDES)

PRIVATE_INCLUDES = pal1.h pal1sofa.h

PUBLIC_C_FILES = \
palAddet.c \
palAirmas.c \
palAltaz.c \
palAmp.c \
palAmpqk.c \
palAop.c \
palAoppa.c \
palAoppat.c \
palAopqk.c \
palAtmdsp.c \
palCaldj.c \
palDafin.c \
palDcmpf.c \
palDe2h.c \
palDeuler.c \
palDfltin.c \
palDh2e.c \
palDjcal.c \
palDmat.c \
palDs2tp.c \
palDat.c \
palDmoon.c \
palDrange.c \
palDt.c \
palDtp2s.c \
palDtps2c.c \
palDtt.c \
palEcleq.c \
palEcmat.c \
palEl2ue.c \
palEpco.c \
palEpv.c \
palEtrms.c \
palEqecl.c \
palEqgal.c \
palEvp.c \
palFitxy.c \
palFk45z.c \
palFk524.c \
palFk54z.c \
palGaleq.c \
palGalsup.c \
palGe50.c \
palGeoc.c \
palIntin.c \
palInvf.c \
palMap.c \
palMappa.c \
palMapqk.c \
palMapqkz.c \
palNut.c \
palNutc.c \
palOap.c \
palOapqk.c \
palObs.c \
palPa.c \
palPcd.c \
palPertel.c \
palPertue.c \
palPlanel.c \
palPlanet.c \
palPlante.c \
palPlantu.c \
palPm.c \
palPolmo.c \
palPrebn.c \
palPrec.c \
palPreces.c \
palPrenut.c \
palPv2el.c \
palPv2ue.c \
palPvobs.c \
palPxy.c \
palRanorm.c \
palRdplan.c \
palRefco.c \
palRefro.c \
palRefv.c \
palRefz.c \
palRverot.c \
palRvgalc.c \
palRvlg.c \
palRvlsrd.c \
palRvlsrk.c \
palSubet.c \
palSupgal.c \
palUe2el.c \
palUe2pv.c \
palUnpcd.c \
palVers.c \
palXy2xy.c \
palZd.c

PRIVATE_C_FILES = \
pal1Atms.c \
pal1Atmt.c

stardocs_DATA = @STAR_LATEX_DOCUMENTATION@

TESTS = palTest

check_PROGRAMS = palTest
palTest_SOURCES = palTest.c
palTest_LDADD = libpal.la

# A target for making the SUN documentation. We do not do this automatically
palsun.tex: $(PUBLIC_C_FILES) palOne2One.c
	-rm -f palsun.tex all.c
	cat $^ > all.c
	${STARCONF_DEFAULT_PREFIX}/bin/sst/prolat in=all.c out=palsun.tex single=no page=no atask=no document=no
	-rm all.c

# Phony target to copy palsun.tex (except palOne2One) into sun267.tex and update the version and date.
updatesun: sun267.tex palsun.tex
	awk '/^%% include palsun\.tex/ {print; print ""; exit}; /^\\stardocversion/ {gsub(/{[0-9.]*}/, "{@PACKAGE_VERSION@}")}; /^\\stardocdate/ {gsub(/{[A-Za-z0-9 ]*}/, strftime("{%Y %B %e}"))}; {print}' $< > sun267_updated.tex
	awk 'BEGIN {onetoone=0}; /^   palOne2One/ {onetoone=1}; (! onetoone) {print}; /^\\sstroutine{/ {onetoone=0}'  palsun.tex >> sun267_updated.tex
	awk 'BEGIN {footer=0}; /^% \? End of main text/{footer=1; print ""; print ""}; footer {print}' $< >> sun267_updated.tex
	mv sun267_updated.tex $<
