#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Intersector requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()
    
# Compilation des fortrans ===================================================
from KCore.config import *
if (f77compiler == "None"):
    print "Error: a fortran 77 compiler is required for compiling Intersector."
args = Dist.getForArgs(); opt = ''
for c in xrange(len(args)):
    opt += 'FOPT'+str(c)+'='+args[c]+' '
os.system("make -e FC="+f77compiler+" WDIR=Intersector/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs and libraries ===========================================
libraryDirs = ["build/"+prod, kcoreLibDir]
libraries = ["intersector", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup ======================================================================
import srcs
setup(
    name="Intersector",
    version="2.5",
    description="Mesh-intersection-based services in *Cassiopee*.",
    author="Onera",
    package_dir={"":"."},
    packages=['Intersector'],
    ext_modules=[Extension('Intersector.intersector',
                           sources=["Intersector/intersector.cpp"]+srcs.cpp_srcs,
                           include_dirs=["Intersector"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
                           library_dirs=additionalLibPaths+libraryDirs,
                           libraries=libraries+additionalLibs,
                           extra_compile_args=Dist.getCppArgs(),
                           extra_link_args=Dist.getLinkArgs()
                           )]
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
