# PAL SConstruct file
import os

version = "0.9.2"

def CheckStarlink(context):
    context.Message( "Checking for Starlink environment...")
    if "STARLINK_DIR" in os.environ:
        star_root = os.environ["STARCONF_DEFAULT_PREFIX"]
        star_lib = os.path.join(star_root, "lib")
        star_inc = os.path.join(star_root, "include")
        star_share = os.path.join(star_root, "share")
        star_bin = os.path.join(star_root, "bin")
        context.env.PrependENVPath('PATH', star_bin)
        context.Result(star_root)
        context.env.Replace(PREFIX = star_root)
        context.env.Append(CPPPATH=star_inc, LIBPATH=star_lib)
        return { "root": star_root,
                 "lib": star_lib,
                 "include": star_inc,
                 "bin": star_bin,
                 "share": star_share, }
    context.Result("failed")
    return None

# Allow --prefix to be specified. We don't use it in Starlink mode
# though
AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix')

env = Environment(PREFIX = "/usr/local" ) # Initialise the environment

prefix = GetOption("prefix")
if prefix:
    env.Replace(PREFIX = prefix)


# Basic configure checks: but not if --help or --clean
if not GetOption("help") and not GetOption("clean"):
    conf = Configure(env, custom_tests = {"CheckStarlink": CheckStarlink} )

    if not conf.CheckCC():
        print("!! Your compiler and/or environment is not correctly configured.")
        Exit(0)

    if conf.CheckFunc("copysign"):
        conf.env.Append(CPPDEFINES={"HAVE_COPYSIGN": 1})

    if conf.CheckFunc("isblank"):
        conf.env.Append(CPPDEFINES={"HAVE_ISBLANK": 1})

    if conf.CheckFunc("strlcpy"):
        conf.env.Append(CPPDEFINES={"HAVE_STRLCPY": 1})

    # Force -lm if we need it
    conf.CheckLib("m","sin")

    # Try to look in current directory
    conf.env.Append(LIBPATH=["."])
    conf.env.Append(CPPPATH=["."])

    # If we are in a Starlink environment we know we have
    # ERFA so just set things up for that. This should be done
    # using a general SCons "are we starlink" plugin
    starlink = conf.CheckStarlink()
    if starlink is not None:
        conf.env.Append(LIBS=["erfa"])
        conf.env.Append(CPPDEFINES={"HAVE_STAR_UTIL": 1})
        conf.env.Append(LIBS=["starutil"])
    else:
        # Allow PREFIX to work
        conf.env.Append(CPPPATH=[os.path.join("$PREFIX", "include")])
        conf.env.Append(LIBPATH=[os.path.join("$PREFIX", "lib")])

        # Maybe starutil will be available
        if conf.CheckLib("starutil"):
            conf.env.Append(CPPDEFINES={"HAVE_STAR_UTIL": 1})

        # Need to look for ERFA vs SOFA
        if not conf.CheckLib("erfa","eraCal2jd"):
            if conf.CheckLib("sofa_c","iauCal2jd"):
                conf.env.Append(CPPDEFINES={"HAVE_SOFA_H": 1})
            else:
                print("!! Neither ERFA not SOFA library located. Can not continue. !!")
                Exit(0)

    env = conf.Finish()

# PAL source code
libpal_sources = [
    "pal1Atms.c",
    "pal1Atmt.c",
    "palAddet.c",
    "palAirmas.c",
    "palAltaz.c",
    "palAmp.c",
    "palAmpqk.c",
    "palAop.c",
    "palAoppa.c",
    "palAoppat.c",
    "palAopqk.c",
    "palAtmdsp.c",
    "palCaldj.c",
    "palDafin.c",
    "palDat.c",
    "palDe2h.c",
    "palDeuler.c",
    "palDfltin.c",
    "palDh2e.c",
    "palDjcal.c",
    "palDmat.c",
    "palDmoon.c",
    "palDrange.c",
    "palDs2tp.c",
    "palDt.c",
    "palDtp2s.c",
    "palDtps2c.c",
    "palDtt.c",
    "palEcmat.c",
    "palEl2ue.c",
    "palEpco.c",
    "palEpv.c",
    "palEqecl.c",
    "palEqgal.c",
    "palEtrms.c",
    "palEvp.c",
    "palFk45z.c",
    "palFk524.c",
    "palFk54z.c",
    "palGaleq.c",
    "palGalsup.c",
    "palGe50.c",
    "palGeoc.c",
    "palIntin.c",
    "palMap.c",
    "palMappa.c",
    "palMapqk.c",
    "palMapqkz.c",
    "palNut.c",
    "palNutc.c",
    "palOap.c",
    "palOapqk.c",
    "palObs.c",
    "palOne2One.c",
    "palPa.c",
    "palPertel.c",
    "palPertue.c",
    "palPlanel.c",
    "palPlanet.c",
    "palPlante.c",
    "palPlantu.c",
    "palPm.c",
    "palPolmo.c",
    "palPrebn.c",
    "palPrec.c",
    "palPreces.c",
    "palPrenut.c",
    "palPv2el.c",
    "palPv2ue.c",
    "palPvobs.c",
    "palRdplan.c",
    "palRefco.c",
    "palRefro.c",
    "palRefv.c",
    "palRefz.c",
    "palRverot.c",
    "palRvgalc.c",
    "palRvlg.c",
    "palRvlsrd.c",
    "palRvlsrk.c",
    "palSubet.c",
    "palSupgal.c",
    "palUe2el.c",
    "palUe2pv.c",
    ]

sun267_sources = [ "sun267.tex" ]
sun267_pdf = env.PDF( sun267_sources )
Default(sun267_pdf)

# Compiler should look in current directory for header files

staticpal = env.StaticLibrary(target="pal", source=libpal_sources)
sharedpal = env.SharedLibrary( target="pal", source = libpal_sources,
                               SHLIBVERSION=version )

palTest = env.Program("palTest", "palTest.c", LIBS=["pal"] )
test_alias = Alias("test", [palTest], palTest[0].abspath)
AlwaysBuild(test_alias)

installed_sharedlib = env.InstallVersionedLib(os.path.join("$PREFIX","lib"),
                                            [sharedpal], SHLIBVERSION=version )
installed_staticlib = env.Install("$PREFIX/lib", [staticpal] )

# Just build the library by default
Default(sharedpal)
Default(staticpal)

# install on request
if "install" in COMMAND_LINE_TARGETS:
    env.Alias( "install", "$PREFIX" )
    Default(installed_sharedlib, installed_staticlib)

