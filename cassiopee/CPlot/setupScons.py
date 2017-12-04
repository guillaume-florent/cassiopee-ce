#!/usr/bin/env python
import os, sys
from distutils.core import setup, Extension
import KCore.config

#=============================================================================
# CPlot requires:
# C++ compiler
# Numpy
# KCore
# GL
# optional: GLEW, PNG, MPEG, OSMesa
#=============================================================================

# If you want to use CPlot as a offscreen plotter (as on clusters)
# set UseOSMesa to True (requires mesa)
UseOSMesa = KCore.config.CPlotOffScreen

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

mySystem = Dist.getSystem()
if (mySystem[0] == 'mingw' and mySystem[1] == '32'):
    libraries = ["cplot", "kcore", "wsock32", "winmm", "gdi32"]
    libGL = ['opengl32', 'glu32']
elif (mySystem[0] == 'mingw' and mySystem[1] == '64'):
    libraries = ["cplot", "kcore", "wsock32", "winmm", "gdi32"]
    libGL = ['opengl32', 'glu32']
elif (mySystem[0] == 'Darwin'):
    libraries = ["kcore", "X11", "Xmu", "cplot"]
    libGL = ['GL', 'GLU'] 
else:
    libraries = ["cplot", "kcore", "Xi", "Xmu", "rt"]
    libGL = ['GL', 'GLU']

from KCore.config import *
if not UseOSMesa: libraries += libGL

prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
libraryDirs = ["build/"+prod]
includeDirs = [numpyIncDir, kcoreIncDir]

(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Test if PNG exists =========================================================
(png, pngIncDir, pngLib) = Dist.checkPng(additionalLibPaths,
                                         additionalIncludePaths)
if png:
    libraries += ["png"]
    if mySystem[0] == 'mingw':
        if Dist.useStatic() == False: libraries += ["zlib1"]
        else: libraries += ["z"]
    libraryDirs += [pngLib]
    includeDirs += [pngIncDir]

# Test if MPEG exists =========================================================
(mpeg, mpegIncDir, mpegLib) = Dist.checkMpeg(additionalLibPaths,
                                             additionalIncludePaths)
if mpeg:
    libraries += ["avcodec"]
    libraryDirs += [mpegLib]
    includeDirs += [mpegIncDir]

# Test if OSMesa exists =======================================================
# Put this to True for using CPlot in batch mode
if UseOSMesa:
    (OSMesa, OSMesaIncDir, OSMesaLib) = Dist.checkOSMesa(additionalLibPaths,
                                                         additionalIncludePaths)
    if OSMesa:
        libraries += ["OSMesa", "GL", "GLU"]
        libraryDirs += [OSMesaLib]
        includeDirs += [OSMesaIncDir]
else: OSMesa = False
    
libraryDirs += [kcoreLibDir]

# Extensions =================================================================
import KCore.installPath
EXTRA = ['-D__SHADERS__']

if OSMesa: EXTRA += ['-D__MESA__']

EXTRA += Dist.getCppArgs()

extensions = [
    Extension('CPlot.cplot',
              sources=['CPlot/cplot.cpp'],
              include_dirs=["CPlot"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=EXTRA,
              extra_link_args=Dist.getLinkArgs()
	)
    ]

# Setup ======================================================================
setup(
    name="CPlot",
    version="2.5",
    description="A plotter for *Cassiopee* Modules.",
    author="Onera",
    package_dir={"":"."},
    packages=['CPlot'],
    ext_modules=extensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
