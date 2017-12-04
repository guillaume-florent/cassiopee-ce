# Functions used in *Cassiopee* modules setup.py
import os, sys, distutils.sysconfig, platform, glob, subprocess

# Toggle to True for compiling for debug (valgrind, inspector)
DEBUG = False
    
#==============================================================================
# Check module import
# Write SUCCESS or FAILED (with colored output)
#==============================================================================
def checkModuleImport(moduleName):
    # Remove . from PYTHONPATH
    try:
        del sys.path[sys.path.index('')]
    except: pass
    try:
        del sys.path[sys.path.index(os.getcwd())]
    except: pass
    # Try to detect if colored output is supported
    color = False
    if sys.stdout.isatty(): color = True

    # try import module
    try:
        __import__(moduleName)
        if color:
            print "\033[32m%s correctly installed.\033[0m"%moduleName
        else: print "%s correctly installed."%moduleName
    except Exception, inst:
        if color:
            print "\033[31mFAILED: %s\033[0m"%inst
            print "\033[31mFAILED: %s badly installed.\033[0m"%moduleName
        else:
            print "FAILED: %s"%inst
            print "FAILED: %s badly installed."%moduleName

#==============================================================================
# Return informations on the current operating system
# Return: Unix, Windows, Darwin, Java, mingw + bits of the system ('32', '64')
#==============================================================================
def getSystem():
     # Detection Mingw
     try:
          # On se base en premier sur ELSAPROD
          key = os.environ['ELSAPROD']
          if key == 'win64': return ['mingw', '64']
          elif key == 'win32': return ['mingw', '32']
          key = os.environ['MSYSTEM']
          if key == 'MINGW32': return ['mingw', '32']
          elif key == 'MINGW64': return ['mingw', '64']
     except: pass
     # System bits
     bits = '32'
     try:
          bits = platform.architecture()[0]
          bits = bits[0:2]
     except: pass
     # Windows, unix, mac
     name = platform.uname()
     return [name[0], bits]

#==============================================================================
# Get name in environ
# Return '' if name is not in environ
#==============================================================================
def getenv(name):
     if os.environ.has_key(name): return os.environ[name]
     else: return ''

#==============================================================================
# Check all
#==============================================================================
def checkAll():
    out = []
    try:
        (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = checkPython()
        out += ['Python: %s'%pythonVersion]
    except: out += ['Python: include is missing.']
    try:
        (numpyVersion, numpyIncDir, ti) = checkNumpy()
        out += ['Numpy: %s'%numpyVersion]
    except: out += ['Numpy: is missing or numpy includes are missing.']
    (ok, CppLibs, CppLibPaths) = checkCppLibs()
    if ok: out += ['C++: OK, %s %s'%(CppLibs, CppLibPaths)]
    else: out += ['C++: Fail']
    (ok, FLibs, FLibPaths) = checkFortranLibs()
    if ok: out += ['f77: OK, %s %s'%(FLibs, FLibPaths)]
    else: out += ['f77: Fail']
    (ok, hdfIncDir, hdfLib) = checkHdf()
    if ok: out += ['hdf: OK, %s %s'%(hdfIncDir, hdfLib)]
    else: out += ['hdf: missing']
    for i in out: print i

#==============================================================================
# Check python includes / libs
#==============================================================================
def checkPython():
    pythonVersion = distutils.sysconfig.get_python_version()
    vars = distutils.sysconfig.get_config_vars()
    
    pythonIncDir = distutils.sysconfig.get_python_inc()
    if not os.path.exists(pythonIncDir):
         raise SystemError("Error: Python includes are required for the compilation of Cassiopee modules.")
        
    pythonLibDir = distutils.sysconfig.get_python_lib()
    try: 
         a = distutils.sysconfig.get_config_var('LDLIBRARY')
         a = a.replace('lib', '')
         a = a.replace('.a', '')
         a = a.replace('.so', '')
         pythonLibs = [a]
    except: pythonLibs = []
    return (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs)

#=============================================================================
# Check numpy
#=============================================================================
def checkNumpy():
     numpyVersion = False
     numpyIncDir = ''
     try:
          import numpy
          numpyIncDir = numpy.get_include()
          numpyVersion = numpy.__version__
     except ImportError:
          raise SystemError("Error: numpy is required for the compilation of Cassiopee modules.")
     
     if not os.path.exists(numpyIncDir):
          raise SystemError("Error: numpy includes are required for the compilation of Cassiopee modules.")
        
     return (numpyVersion, numpyIncDir, '')

#=============================================================================
# Retourne le chemin d'installation des modules comme cree par distUtils
#=============================================================================
def getInstallPath(prefix):
    mySystem = getSystem()[0]; bits = getSystem()[1]
    # Based on spec
    if mySystem == 'Windows' or mySystem == 'mingw':
         installPath = prefix + "/Lib/site-packages"
    elif mySystem == 'Darwin':
         pythonLib = distutils.sysconfig.get_python_lib()
         pythonLib = pythonLib.split('/')
         pythonVersion = pythonLib[-2]
         installPath = prefix + '/lib/python'+pythonVersion+'/site-packages'
    else: # unix
         pythonLib = distutils.sysconfig.get_python_lib()
         pythonLib = pythonLib.split('/')
         # Based on python lib
         #installPath = prefix + '/' + '/'.join(pythonLib[-3:])
         # Python version
         pythonVersion = pythonLib[-2]
         Site = pythonLib[-1]
         # Lib
         Lib = pythonLib[-3]
         #if (bits == '64'): Lib = 'lib64'
         installPath = '%s/%s/%s/site-packages'%(prefix, Lib, pythonVersion)
    return installPath

#=============================================================================
# Write installPath, the installation path of Cassiopee to installPath.py
#=============================================================================
def writeInstallPath():
    import re
    prefix = sys.prefix
    a = sys.argv
    for i in a:
        if re.compile('--prefix=').search(i) is not None: prefix = i[9:]
        elif re.compile('prefix=').search(i) is not None: prefix = i[7:]
    installPath = getInstallPath(prefix)
    p = open('installPath.py', 'w')
    if p is None:
        raise SystemError("Error: can not open file installPath.py for writing.")
    p.write('installPath = \'%s\'\n'%installPath)
    mySystem = getSystem()[0]; bits = getSystem()[1]
    if mySystem == 'Windows' or mySystem == 'mingw': Lib = 'Lib'
    elif mySystem == 'Darwin': Lib = 'lib'
    else:
         pythonLib = distutils.sysconfig.get_python_lib()
         pythonLib = pythonLib.split('/')
         Lib = pythonLib[-3]
    p.write('libPath = \'%s/%s\'\n'%(prefix,Lib))
    p.write('includePath = \'%s\'\n'%(os.getcwd()))
    p.close()

#==============================================================================
# Write env files
# Directement dans le repertoire d'installation
#==============================================================================
def writeEnvs():
     try: import KCore.installPath as K
     except: import installPath as K
     libPath = K.libPath
     installPathLocal = K.installPath
     env = os.environ.data
     cassiopee = env.get('CASSIOPEE', '')
     elsaprod = env.get('ELSAPROD', '')
     if cassiopee != '': envPath = libPath+'/../../../'
     else: envPath = libPath+'/../'
     cmdPath = libPath+'/..'
     installLD = os.getenv('LD_LIBRARY_PATH')

     # max cores
     try:
          import multiprocessing
          mt = multiprocessing.cpu_count()
     except: mt = 1

     # sh
     p = open(envPath+"env_Cassiopee.sh", 'w')
     if cassiopee != '': p.write("export CASSIOPEE=%s\n"%cassiopee)
     if elsaprod != '': p.write("export ELSAPROD=%s\n"%elsaprod)
     p.write("export OMP_NUM_THREADS=%d\n"%mt)
     p.write("export PATH=%s:$PATH\n"%cmdPath)
     p.write("if [ \"$PYTHONPATH\" = \"\" ]; then\n")
     p.write("      export PYTHONPATH=%s\n"%installPathLocal)
     p.write("else\n")
     p.write("      export PYTHONPATH=%s:$PYTHONPATH\n"%installPathLocal)
     p.write("fi\n")
     if installLD is None:
          p.write("if [ \"$LD_LIBRARY_PATH\" = \"\" ]; then\n")
          p.write("      export LD_LIBRARY_PATH=%s\n"%libPath)
          p.write("else\n")
          p.write("      export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n"%libPath)
          p.write("fi\n")
     else:
          p.write("if [ \"$LD_LIBRARY_PATH\" = \"\" ]; then\n")
          p.write("      export LD_LIBRARY_PATH=%s:%s\n"%(libPath,installLD))
          p.write("else\n")
          p.write("      export LD_LIBRARY_PATH=%s:%s:$LD_LIBRARY_PATH\n"%(libPath,installLD))
          p.write("fi\n")
     p.close()
     
     # csh
     p = open(envPath+"env_Cassiopee.csh", 'w')
     if cassiopee != '': p.write("setenv CASSIOPEE %s\n"%cassiopee)
     if elsaprod != '': p.write("setenv ELSAPROD %s\n"%elsaprod)
     p.write("setenv OMP_NUM_THREADS %d\n"%mt)
     p.write("set path=(%s $path)\n"%cmdPath)
     p.write("if ($?PYTHONPATH == 0) then\n")
     p.write("     setenv PYTHONPATH %s\n"%installPathLocal)
     p.write("else\n")
     p.write("     setenv PYTHONPATH %s:$PYTHONPATH\n"%installPathLocal)
     p.write("endif\n")
     if installLD is None:
          p.write("if ($?LD_LIBRARY_PATH == 0) then\n")
          p.write("     setenv LD_LIBRARY_PATH %s\n"%libPath)
          p.write("else\n")
          p.write("     setenv LD_LIBRARY_PATH %s:$LD_LIBRARY_PATH\n"%libPath)
          p.write("endif\n")
     else:
          p.write("if ($?LD_LIBRARY_PATH == 0) then\n")
          p.write("     setenv LD_LIBRARY_PATH %s:%s\n"%(libPath,installLD))
          p.write("else\n")
          p.write("     setenv LD_LIBRARY_PATH %s:%s:$LD_LIBRARY_PATH\n"%(libPath,installLD))
          p.write("endif\n")
     p.close()
     
     # bat
     p = open(envPath+"env_Cassiopee.bat", 'w')
     p.write("path = "+libPath+";"+cmdPath+"%PATH%\n")
     p.write("set PYTHONPATH="+installPathLocal+";%PYTHONPATH%\n")
     p.write("set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%\n")
     p.close()
     
#==============================================================================
# Write setup.cfg en fonction du compilateur C++ (si different de None)
# setup.cfg est utilise par setup de python pour choisir le compilo.
#==============================================================================
def writeSetupCfg():
    try: from config import Cppcompiler
    except: from KCore.config import Cppcompiler
    mySystem = getSystem()
    
    # Windows + mingw
    if mySystem[0] == 'mingw' and mySystem[1] == '32':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=mingw32\n')
        p.close(); return
    if mySystem[0] == 'mingw' and mySystem[1] == '64':
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=mingw32\n')
        p.close(); return
    
    # Unix
    if Cppcompiler == "None" or Cppcompiler == "":
        a = os.access("./setup.cfg", os.F_OK)
        if a: os.remove("./setup.cfg")
    elif Cppcompiler.find('icc') == 0 or Cppcompiler.find('icpc') == 0:
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=intel\n')
        p.close()
    elif Cppcompiler.find('gcc') == 0 or Cppcompiler.find('g++') == 0:
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=unix\n')
        p.close()
    else:
        p = open("./setup.cfg", 'w')
        p.write('[build_ext]\ncompiler=%s\n'%Cppcompiler)
        p.close()
        
    if Cppcompiler.find('icc') == 0 or Cppcompiler.find('icpc') == 0:
         import numpy.distutils
         numpy.distutils.ccompiler.new_compiler(compiler='intel')

#==============================================================================
# Retourne le compilo c, c++ et ses options tels que definis dans distutils
# ou dans config.py (installBase.py)
#==============================================================================
def getDistUtilsCompilers():
    vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT',
                                               'BASECFLAGS', 'CCSHARED',
                                               'LDSHARED', 'SO')
    for i in xrange(len(vars)):
        if vars[i] is None: vars[i] = ""

    try: from config import Cppcompiler
    except: from KCore.config import Cppcompiler
    if Cppcompiler != 'None' or Cppcompiler != '':
        if Cppcompiler.find('g++') != -1: # g++-version + mingw-g++-version
            vars[0] = Cppcompiler.replace('g++', 'gcc'); vars[1] = Cppcompiler
        elif Cppcompiler.find('gcc') != -1:
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('gcc', 'g++')
        elif Cppcompiler.find('icpc') != -1:
            vars[0] = Cppcompiler.replace('icpc', 'icc'); vars[1] = Cppcompiler
        elif Cppcompiler.find('icc') != -1:
            vars[0] = Cppcompiler; vars[1] = Cppcompiler.replace('icc', 'icpc')

    (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext) = vars
    cc = cc.split(' ') # enleve les options si mises dans cc
    if len(cc) > 1: cc = cc[0]
    cxx = cxx.split(' ')
    if len(cxx) > 1: cxx = cxx[0]
    return (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext)

#==============================================================================
# Retourne le pre-processeur utilise pour les fichiers fortrans
# IN: config.Cppcompiler
#==============================================================================
def getPP():
    try: from config import Cppcompiler
    except: from KCore.config import Cppcompiler
    sizes = '-DINTEGER_E="INTEGER*4" -DREAL_E="REAL*8"'
    if Cppcompiler == 'icl.exe': PP = 'fpp.exe '+sizes+' \I'
    elif Cppcompiler == "x86_64-w64-mingw32-gcc":
         PP = 'x86_64-w64-mingw32-cpp -traditional %s -I'%sizes
    else: PP = 'cpp -traditional %s -I'%sizes
    return PP

#==============================================================================
# Retourne l'achiveur pour faire les librairies statiques
# IN: config.Cppcompiler
#==============================================================================
def getAR():
    try: from config import Cppcompiler
    except: from KCore.config import Cppcompiler
    if Cppcompiler == "icl.exe": return 'ar.exe '
    elif Cppcompiler == "x86_64-w64-mingw32-gcc":
         return 'x86_64-w64-mingw32-ar'
    else: return 'ar'

#==============================================================================
# Retourne le prefix pour le repertoire ou on stocke les modules f90
# IN: config.f90compiler
#==============================================================================
def getFortranModDirPrefix():
    try: from config import f90compiler
    except: from KCore.config import f90compiler
    if f90compiler == 'ifort': return '-module'
    elif f90compiler == 'ifort.exe': return '-module'
    elif f90compiler == 'gfortran': return '-J'
    elif f90compiler == 'g95': return '-fmod'
    else: return ''

#==============================================================================
# Retourne 1 si oui
# IN: config.useOMP
#==============================================================================
def useOMP():
     try: from config import useOMP
     except: from KCore.config import useOMP
     if useOMP: return 1
     else: return 0

#==============================================================================
# Retourne 1 si on produit des librairies statiques
# IN: config.useStatic
#==============================================================================
def useStatic():
     try: from config import useStatic
     except: from KCore.config import useStatic
     if useStatic: return 1
     else: return 0

#==============================================================================
# Retourne les versions des compilateurs
#==============================================================================
def getVersion(compiler):
     major = 0; minor = 0
     cmd = [compiler, '--version']
     proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
     out = ''
     while True:
          line = proc.stdout.readline()
          if line != '': out += line
          else: break
     out = out.split('\n')
     out = out[0]
     out = out.split(' ')
     for i in out:
          isVersion = i.split('.')
          if len(isVersion)>1: # maybe
               try: 
                    major = int(isVersion[0])
                    minor = int(isVersion[1])
                    break
               except: continue
     return (major, minor)

def getCppVersion():
     try: from config import Cppcompiler
     except: from KCore.config import Cppcompiler
     return getVersion(Cppcompiler)

def getForVersion():
     try: from config import f77compiler
     except: from KCore.config import f77compiler
     return getVersion(f77compiler)

#=============================================================================
# Retourne le nbre de lignes du cache
# Se base sur l'option du compilateur C si elle contient -DCACHELINE=XX
#==============================================================================
def getCacheLine():
     opts = getCppArgs()
     for i in opts:
          if i[0:11] == '-DCACHELINE':
               val = int(i[12:]); return val
     return 1

#==============================================================================
# Retourne les options additionelles des compilos definies dans config
#==============================================================================
def getCppAdditionalOptions():
     try: from config import CppAdditionalOptions
     except: from KCore.config import CppAdditionalOptions
     return CppAdditionalOptions

def getf77AdditionalOptions():
     try: from config import f77AdditionalOptions
     except: from KCore.config import f77AdditionalOptions
     return f77AdditionalOptions

#==============================================================================
# Retourne les arguments pour le compilateur Cpp
# IN: config.Cppcompiler, config.useStatic, config.useOMP, 
# config.CppAdditionalOptions
#==============================================================================
def getCppArgs():
    try: from config import Cppcompiler
    except: from KCore.config import Cppcompiler
    mySystem = getSystem()
    compiler = Cppcompiler.split('/')
    l = len(compiler)-1
    Cppcompiler = compiler[l]
    if Cppcompiler == "None": return []
    options = getCppAdditionalOptions()
    if Cppcompiler.find("icpc") == 0 or Cppcompiler.find("icc") == 0: 
         if DEBUG: options += ['-g', '-O0', '-wd47', '-wd1224']
         else: options += ['-DNDEBUG', '-O2', '-wd47', '-wd1224']
         v = getCppVersion()
         if v[0] < 15: options += ['-fp-speculation=strict']
         else: options += ['-fp-model=strict']
         if useOMP() == 1:
            if v[0] < 15: options += ['-openmp']
            else: options += ['-qopenmp']
         if useStatic() == 1: options += ['-static']
         else: options += ['-fPIC']
         return options
    elif Cppcompiler.find("gcc") == 0 or Cppcompiler.find("g++") == 0:
         if DEBUG: options += ['-g', '-O0', '-Wall']
         else: options += ['-DNDEBUG', '-O3', '-Wall']
         if useOMP() == 1: options += ['-fopenmp']
         if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
         else: options += ['-fPIC']
         if mySystem[0] == 'mingw' and mySystem[1] == '32':
              options.remove('-fPIC')
              options += ['-large-address-aware']
         return options
    elif Cppcompiler == "icl.exe":
         options += ['/EHsc', '/MT']
         if useOMP() == 1: options += ['/Qopenmp']
         return options
    elif Cppcompiler == "x86_64-w64-mingw32-gcc" or Cppcompiler == "x86_64-w64-mingw32-g++":
         options += ['-DMS_WIN64', '-fpermissive', '-D__USE_MINGW_ANSI_STDIO=1']
         if DEBUG: options += ['-g', 'O0']
         else: options += ['-DNDEBUG', '-O3']
         if useOMP() == 1: options += ['-fopenmp']
         if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
         else: options += ['-fPIC']
         return options
    elif Cppcompiler.find("clang") == 0 or Cppcompiler.find("clang++") == 0:
         if DEBUG: options += ['-g', '-O0', '-Wall']
         else: options += ['-DNDEBUG', '-O3', '-Wall']
         if useOMP() == 1: options += ['-fopenmp']
         if useStatic() == 1: options += ['--static', '-static-libstdc++', '-static-libgcc']
         else: options += ['-fPIC']
         return options
    else: return options

#==============================================================================
# Retourne les arguments pour le compilateur Fortran
# IN: config.f77compiler
#==============================================================================
def getForArgs():
    try: from config import f77compiler
    except: from KCore.config import f77compiler
    mySystem = getSystem()
    compiler = f77compiler.split('/')
    l = len(compiler)-1
    f77compiler = compiler[l]
    if f77compiler == "None": return []
    options = getf77AdditionalOptions()
    if f77compiler.find("gfortran") == 0:
         if DEBUG: options += ['-g', '-O0']
         else: options += ['-O3']
         if useOMP() == 1: options += ['-fopenmp']
         if useStatic() == 1: options += ['--static']
         else: options += ['-fPIC']
         if (mySystem[0] == 'mingw' and mySystem[1] == '32'):
              options.remove('-fPIC')
              options += ['-large-address-aware']
         return options
    elif f77compiler.find("ifort") == 0:
         if DEBUG: options += ['-g', '-O0', '-fPIC']
         else: options += ['-O3']
         if useOMP() == 1: 
            v = getForVersion()
            if v[0] < 15: options += ['-openmp']
            else: options += ['-qopenmp']
         if useStatic() == 1: options += ['-static']
         else: options += ['-fPIC']
         return options
    elif f77compiler == "x86_64-w64-mingw32-gfortran":
         if DEBUG: options += ['-g', '-O0']
         else: options += ['-O3']
         if useOMP() == 1: options += ['-fopenmp']
         if useStatic() == 1: options += ['--static']
         else: options += ['-fPIC']
         return options
    elif f77compiler == "ifort.exe":
         if useOMP() == 1: return ['/names:lowercase', '/assume:underscore', '/Qopenmp']
         else: return ['/names:lowercase', '/assume:underscore']
    else: return options

#==============================================================================
# Retourne les arguments pour le linker
# IN: config.Cppcompiler, config.useStatic, config.useOMP
#==============================================================================
def getLinkArgs():
     try: from config import Cppcompiler
     except: from KCore.config import Cppcompiler
     out = []
     if Cppcompiler == 'gcc' or Cppcompiler == 'g++':
          if useStatic() == 1: out += ['--static']
     elif Cppcompiler == 'icc':
          if useStatic() == 1: out += ['-static']
     elif Cppcompiler == "x86_64-w64-mingw32-gcc":
          if useStatic() == 1: out += ['--static']
     mySystem = getSystem()[0]
     if mySystem == 'Darwin':
	  if useStatic() == 0: out += ['-dynamiclib'] 
     return out

#=============================================================================
# Check PYTHONPATH
# Verifie que installPath est dans PYTHONPATH
#=============================================================================
def checkPythonPath():
    import re
    try: import KCore.installPath as K
    except: import installPath as K
    installPathLocal = K.installPath
    a = os.getenv("PYTHONPATH")
    if a is None:
        print 'Warning: to use the module, please add: '\
        + installPathLocal + ' to your PYTHONPATH.'
    else:
        if (re.compile(installPathLocal).search(a) is None):
            print 'Warning: to use the module, please add: '\
            + installPathLocal + ' to your PYTHONPATH.'

#=============================================================================
# Check LD_LIBRARY_PATH
#=============================================================================
def checkLdLibraryPath():
    import re
    try: import KCore.installPath as K
    except: import installPath as K
    libPath = K.libPath
    a = os.getenv("LD_LIBRARY_PATH")
    b = os.getenv("LIBRARY_PATH")
    if a is None and b is None:
        print "Warning: to use the module, please add: %s to your LD_LIBRARY_PATH (unix) or PATH (windows)."%libPath
    else:
         if a is not None: ret = a
         else: ret = b
         if (re.compile(libPath).search(ret) is None):
              print "Warning: to use the module, please add: %s to your LD_LIBRARY_PATH (unix) or PATH (windows)."%libPath

#=============================================================================
# Check for KCore
#=============================================================================
def checkKCore():
    try:
        import KCore
        import KCore.installPath
        kcoreIncDir = KCore.installPath.includePath
        kcoreIncDir = os.path.join(kcoreIncDir, 'KCore')
        kcoreLib = KCore.installPath.libPath
        return (KCore.__version__, kcoreIncDir, kcoreLib)

    except ImportError:
        raise SystemError("Error: kcore library is required for the compilation of this module.")
        
#=============================================================================
# Check for Generator
#=============================================================================
def checkGenerator():
    try:
        import Generator
        import KCore.installPath
        generatorIncDir = KCore.installPath.includePath
        generatorIncDir = os.path.dirname(generatorIncDir)
        generatorIncDir = os.path.join(generatorIncDir, 'Generator/Generator')
        generatorLib = KCore.installPath.libPath
        return (Generator.__version__, generatorIncDir, generatorLib)

    except ImportError:
        raise SystemError("Error: generator library is required for the compilation of this module.")

#=============================================================================
# Check for Fast module
#=============================================================================
def checkFast():
    try:
        import Fast
        import KCore.installPath
        fastIncDir = KCore.installPath.includePath
        fastIncDir = os.path.dirname(fastIncDir)
        fastIncDir = os.path.join(fastIncDir, 'Fast/Fast')
        fastLib = KCore.installPath.libPath
        return (Fast.__version__, fastIncDir, fastLib)

    except ImportError:
        raise SystemError("Error: fast library is required for the compilation of this module.")

#=============================================================================
# Check for Cassiopee Kernel in Dist/ or Kernel/
# La prod de Cassiopee ne doit pas contenir mpi
#=============================================================================
def checkCassiopee():
    Cassiopee = False
    CassiopeeLibDir = ""
    CassiopeeIncDir = ""
    CassiopeeUseMpi = False # si le Kernel utilise mpi

    kvar = os.getenv("CASSIOPEE")
    pvar = os.getenv("ELSAPROD") # prod de Cassiopee
    if kvar is None:
        return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)
    if pvar is None:
        return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)

    # Cassiopee must exists and be installed
    # si Kernel existe, on l'utilise en priorite
    # sinon, on utilise Dist
    a1 = os.access(kvar+"/Kernel/include", os.F_OK)
    a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsAc.so', os.F_OK)
    b1 = os.access(kvar+"/Dist/include", os.F_OK)
    b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.so', os.F_OK)

    if (a1 and a2):
        Cassiopee = True
        CassiopeeIncDir = kvar + "/Kernel/include"
        CassiopeeLibDir = kvar + "/Kernel/lib/" + pvar

    elif (b1 and b2):
        Cassiopee = True
        CassiopeeIncDir = kvar + "/Dist/include"
        CassiopeeLibDir = kvar + "/Dist/bin/" + pvar

    if not Cassiopee: # essai avec une production mpi
         pvar = pvar.split('_')
         if len(pvar) >= 2:
              pvar = pvar[0]+'_mpi_'+pvar[1]
              a1 = os.access(kvar+"/Kernel/include", os.F_OK)
              a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsAc.so', os.F_OK)
              b1 = os.access(kvar+"/Dist/include", os.F_OK)
              b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.so', os.F_OK)
              if (a1 and a2):
                   Cassiopee = True
                   CassiopeeUseMpi = True
                   CassiopeeIncDir = kvar + "/Kernel/include"
                   CassiopeeLibDir = kvar + "/Kernel/lib/" + pvar

              elif (b1 and b2):
                   Cassiopee = True
                   CassiopeeUseMpi = True
                   CassiopeeIncDir = kvar + "/Dist/include"
                   CassiopeeLibDir = kvar + "/Dist/bin/" + pvar

    if Cassiopee:
        print 'Info: Cassiopee Kernel detected at '+CassiopeeLibDir+'.'
        print 'Info: .Cassiopee extension will be built.'

    return (Cassiopee, CassiopeeIncDir, CassiopeeLibDir, CassiopeeUseMpi)

#=============================================================================
# Check for elsA Kernel in Dist/ or Kernel/
#=============================================================================
def checkElsa():
    elsA = False
    elsALibDir = ""
    elsAIncDir = ""
    elsAUseMpi = False

    kvar = os.getenv("ELSA")
    if kvar is None: return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)
    pvar = os.getenv("ELSAPROD")
    if pvar is None: return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)

    a1 = os.access(kvar+"/Kernel/include", os.F_OK)
    a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsA.x', os.F_OK)
    b1 = os.access(kvar+"/Dist/include", os.F_OK)
    b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsA.x', os.F_OK)
    if (a1 and a2):
        elsA = True
        elsAIncDir = kvar + "/Kernel/include"
        elsALibDir = kvar + "/Kernel/lib/" + pvar

    elif (b1 and b2):
        elsA = True
        elsAIncDir = kvar + "/Dist/include"
        elsALibDir = kvar + "/Dist/bin/" + pvar
    
    if not elsA: # essai avec une production mpi
         pvar = pvar.split('_')
         if len(pvar) >= 2:
              pvar = pvar[0]+'_mpi_'+pvar[1]
              a1 = os.access(kvar+"/Kernel/include", os.F_OK)
              a2 = os.access(kvar+"/Kernel/lib/"+pvar+'/elsA.x', os.F_OK)
              b1 = os.access(kvar+"/Dist/include", os.F_OK)
              b2 = os.access(kvar+"/Dist/bin/"+pvar+'/elsAc.x', os.F_OK)
              if (a1 and a2):
                   elsA = True
                   elsAUseMpi = True
                   elsAIncDir = kvar + "/Kernel/include"
                   elsALibDir = kvar + "/Kernel/lib/" + pvar

              elif (b1 and b2):
                   elsA = True
                   elsAUseMpi = True
                   elsAIncDir = kvar + "/Dist/include"
                   elsALibDir = kvar + "/Dist/bin/" + pvar

    if elsA:
        print 'Info: elsA Kernel detected at '+elsALibDir+'.'
        print 'Info: .Elsa extension will be built.'
    return (elsA, elsAIncDir, elsALibDir, elsAUseMpi)

#=============================================================================
# Check for Glut (libglut)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
#=============================================================================
def checkGlut(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libglut.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libglut.a', additionalLibPaths)
        if l is None:
             l = checkLibFile__('libfreeglut.a', additionalLibPaths)
    i = checkIncFile__('GL/glut.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: glut detected at '+l+'.'
        return (True, i, l)
    else:
        print 'Info: libglut or GL/glut.h was not found on your system.'
        return (False, '', '')

#=============================================================================
# Check for Glew (libglew)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkGlew(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libGLEW.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libGLEW.a', additionalLibPaths)
        if l is None:
             l = checkLibFile__('libglew32.a', additionalLibPaths)
    i = checkIncFile__('GL/glew.h', additionalIncludePaths)
    
    if (i is not None and l is not None):
        print 'Info: glew detected at '+l+'.'
        return (True, i, l)
    else:
        print 'Info: libglew or GL/glew.h was not found on your system. No shader support for CPlot.'
        return (False, '', '')

#=============================================================================
# Check for osmesa (offline rendering)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkOSMesa(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libOSMesa.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libOSMesa.a', additionalLibPaths)
    i = checkIncFile__('GL/osmesa.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: libOSmesa detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libOSMesa or GL/osmesa.h was not found on your system. No offscreen support for CPlot.'
        return (False, '', '')

#=============================================================================
# Check for png (libpng)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkPng(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libpng.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libpng.a', additionalLibPaths)
        if l is None:
             l = checkLibFile__('libpng.dll.a', additionalLibPaths)
    i = checkIncFile__('png.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: png detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libpng or png.h was not found on your system. No png support.'
        return (False, '', '')

#=============================================================================
# Check for mpeg (ffmpeg)
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkMpeg(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libavcodec.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libavcodec.a', additionalLibPaths)
        if l is None:
             l = checkLibFile__('libavcodec.dll.a', additionalLibPaths)
    i = checkIncFile__('libavcodec/avcodec.h', additionalIncludePaths)
    if i is not None: 
         i = checkIncFile__('libavutil/mem.h', additionalIncludePaths)
    if i is not None: 
         i = checkIncFile__('libavutil/imgutils.h', additionalIncludePaths)     
    if (i is not None and l is not None):
        print 'Info: mpeg detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libavcodec or libavcodec/avcodec.h,  libavutil/mem.h or libavutil/imgutils.h was not found on your system. No mpeg support.'
        return (False, '', '')
    
#=============================================================================
# Check for Adf
# additionalPaths: chemins d'installation non standards: ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkAdf(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libcgns.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libcgns.a', additionalLibPaths)
    i = checkIncFile__('adf/ADF.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: Adf detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libadf or adf/ADF.h was not found on your system. No adf support.'
        return (False, '', '')  

#=============================================================================
# Check for Hdf
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkHdf(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libhdf5.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libhdf5.a', additionalLibPaths)
    i = checkIncFile__('hdf5.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: Hdf5 detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libhdf5 or hdf5.h was not found on your system. No hdf5 support.'
        return (False, '', '')

#=============================================================================
# Check for Mpi
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkMpi(additionalLibPaths=[], additionalIncludePaths=[]):
    l = checkLibFile__('libmpi.so', additionalLibPaths)
    if l is None:
        l = checkLibFile__('libmpi.a', additionalLibPaths)
    i = checkIncFile__('mpi.h', additionalIncludePaths)
    if (i is not None and l is not None):
        print 'Info: Mpi detected at %s.'%l
        return (True, i, l)
    else:
        print 'Info: libmpi or mpi.h was not found on your system. No Mpi support.'
        return (False, '', '')

#=============================================================================
# Check for Mpi4py
# additionalPaths: chemins d'installation non standards : ['/home/toto',...]
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkMpi4py(additionalLibPaths=[], additionalIncludePaths=[]):
    try: import mpi4py
    except: return (False, '', '')

    incPaths = []
    try: import KCore.installPath as K
    except: import installPath as K
    incPaths += [K.installPath+'/mpi4py/include']
    import mpi4py
    fileN = mpi4py.__file__
    incPaths += [os.path.dirname(fileN)+'/include']
    i = checkIncFile__('mpi4py/mpi4py.MPI.h', additionalIncludePaths+incPaths)

    if i is not None:
        print 'Info: Mpi4py detected at %s.'%i
        return (True, i, '')
    else:
        print 'Info: mpi4py or mpi4py.MPI.h was not found on your system. No Mpi support.'
        return (False, '', '')

#=============================================================================
# Check for Cython
# Retourne: (True/False, chemin des includes, chemin de la librairie)
#=============================================================================
def checkCython(additionalLibPaths=[], additionalIncludePaths=[]):
     try: 
          import Cython.Compiler.Main as cython_compiler
          return True
     except:
          return False

def cythonize(src):
     import Cython.Compiler.Main as cython_compiler
     sys.stderr.write("cythonize: %r\n" % (src,))
     cython_compiler.compile([src], cplus=True)

#=============================================================================
# Check fortran libs
# additionalLibs: noms des libraries utilises par fortran si non conventionnels
# additionalLibPaths: chemins d'installation des libraries non standards: ['/home/toto',...]
# Retourne (True, [librairies utiles pour le fortran], [paths des librairies])
#=============================================================================
def checkFortranLibs(additionalLibs=[], additionalLibPaths=[],
                     f77compiler=None, useOMP=None):
     if f77compiler is None:
          try: from config import f77compiler
          except: 
            try: from KCore.config import f77compiler
            except: f77compiler = 'gfortran'
     if useOMP is None:
          try: from config import useOMP
          except: 
            try: from KCore.config import useOMP
            except: useOMP = True
     ret = True; libs = []; paths = []
     
     # librairies speciales (forcees sans check)
     libs += additionalLibs
     paths += additionalLibPaths
               
     # gfortran (gfortran, gomp)
     if f77compiler == 'gfortran':
          l = checkLibFile__('libgfortran.so*', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libgfortran.a', additionalLibPaths)
               
          if l is not None:
               libs += ['gfortran']; paths += [l]

          if useOMP:
               l = checkLibFile__('libgomp.so', additionalLibPaths)
               if l is None:
                    l = checkLibFile__('libgomp.a', additionalLibPaths)
               if l is not None:
                    libs += ['gomp']; paths += [l]
               else: ret = False

     # ifort (ifcore, svml, irc, guide, iomp5) 
     if f77compiler == 'ifort':
          l = checkLibFile__('libifcore.so*', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libifcore.a', additionalLibPaths)
               
          if l is not None:
               libs += ['ifcore']; paths += [l]

          l = checkLibFile__('libsvml.so*', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libsvml.a', additionalLibPaths)
          if l is not None:
               libs += ['svml']; paths += [l]

          l = checkLibFile__('libirc.so*', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libirc.a', additionalLibPaths)
          if l is not None:
               libs += ['irc']; paths += [l]

          if useOMP:
               # check guide
               l = checkLibFile__('libguide.so*', additionalLibPaths)
               if l is None:
                    l = checkLibFile__('libguide.a', additionalLibPaths)
               if l is not None:
                    libs += ['guide']; paths += [l]
               # check iomp5
               if l is None:
                    l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                    if l is None:
                         l = checkLibFile__('libiomp5.a', additionalLibPaths)
                    if l is not None:
                         libs += ['iomp5']; paths += [l]
                    else: ret = False
     return (ret, libs, paths)

#=============================================================================
# Check Cpp libs
# additionalLibs: si les libs requises ont des noms non conventionnels
# additionalLibPaths: si chemins des libraries necessaires non standards : ['/home/toto',...]
# Retourne (True, [librairies utiles a cpp], [paths des librairies])
#=============================================================================
def checkCppLibs(additionalLibs=[], additionalLibPaths=[], Cppcompiler=None,
                 useOMP=None):
     if Cppcompiler is None:
          try: from config import Cppcompiler
          except: 
            try: from KCore.config import Cppcompiler
            except: Cppcompiler = 'gcc'
     if useOMP is None:
          try: from config import useOMP
          except: 
            try: from KCore.config import useOMP
            except: useOMP = True

     ret = True; libs = []; paths = []
     
     # librairies additionales (forcees sans check)
     libs += additionalLibs
     paths += additionalLibPaths
               
     # gcc (stdc++, gomp)
     if Cppcompiler == 'gcc' or Cppcompiler == 'g++':
          l = checkLibFile__('libstdc++.so', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libstdc++.a', additionalLibPaths)
               
          if l is not None:
               libs += ['stdc++']; paths += [l]

          if useOMP:
               l = checkLibFile__('libgomp.so', additionalLibPaths)
               if l is None:
                    l = checkLibFile__('libgomp.a', additionalLibPaths)
               if l is not None:
                    libs += ['gomp']; paths += [l]
               else: ret = False
               
     # icc (stdc++, guide ou iomp5)
     if Cppcompiler == 'icc':
          l = checkLibFile__('libstdc++.so*', additionalLibPaths)
          if l is None:
               l = checkLibFile__('libstdc++.a', additionalLibPaths)
               
          if l is not None:
               libs += ['stdc++']; paths += [l]

          if useOMP:
               l = checkLibFile__('libguide.so*', additionalLibPaths)
               if l is None:
                    l = checkLibFile__('libguide.a', additionalLibPaths)
               if l is not None:
                    libs += ['guide']; paths += [l]
               if l is None:
                    l = checkLibFile__('libiomp5.so*', additionalLibPaths)
                    if l is None:
                         l = checkLibFile__('libiomp5.a', additionalLibPaths)
                    if l is not None:
                         libs += ['iomp5']; paths += [l]
                    else: ret = False
     return (ret, libs, paths)

#==============================================================================
# Check if file exists in std lib dirs, path dirs and specified lib dirs
# additionalPaths est en premier pour pouvoir forcer une librairie par config.
#==============================================================================
def checkLibFile__(file, additionalLibPaths):
    p = []
    for i in additionalLibPaths: p.append(i)
    p += ['/usr/local/lib', '/opt/lib', '/usr/lib', '/opt/local/lib']
    mySystem = getSystem()
    env = os.environ
    if mySystem[0] == 'Windows':
         p1 = env.get('PATH', None)
         if p1 is not None: p += p1.split(';')
    else: # unix, mingw...
         p1 = env.get('PATH', None)
         if p1 is not None: p += p1.split(';')
         p1 = env.get('LD_LIBRARY_PATH', None)
         if p1 is not None: p += p1.split(':')
    for i in p:
         a = glob.glob(i+'64/'+file)
         if a != []: return i+'64'
         a = glob.glob(i+'/'+file)
         if a != []: return i
    return None

#==============================================================================
# Check if file exists in std include dir, path dirs and specified include dirs
#==============================================================================
def checkIncFile__(file, additionalIncludePaths):
    p = []
    for i in additionalIncludePaths: p.append(i)
    p += ['/usr/local/include', '/opt/include', '/usr/include', 
          '/opt/local/include']
    mySystem = getSystem()
    env = os.environ
    pp = []
    if mySystem[0] == 'Windows':
         p1 = env.get('PATH', None)
         if p1 is not None: pp += p1.split(';')
    else: # unix, mingw...
         p1 = env.get('PATH', None)
         if p1 is not None: pp += p1.split(':')
         p1 = env.get('LD_LIBRARY_PATH', None)
         if p1 is not None: pp += p1.split(':')
    for i in xrange(len(pp)): pp[i] = pp[i].replace('lib', 'include')
    p += pp
    
    for i in p:
        a = os.access(i+'/'+file, os.F_OK)
        if a: return i
    return None

#==============================================================================
# Ecrit les infos de build dans un fichier buildInfo.py
#==============================================================================
def writeBuildInfo():
     p = open("buildInfo.py", 'w')
     if p is None:
        raise SystemError("Error: can not open file buildInfo.py for writing.")
     try: import KCore.config as config
     except: import config

     dict = {}
     # Date
     import time, pickle
     execTime = time.strftime('%d/%m/%y %Hh%M', time.localtime())
     dict['date'] = execTime

     # Check python
     (pythonVersion, pythonIncDir, pythonLibDir, pythonLibs) = checkPython()
     if (pythonVersion != False): dict['python'] = pythonVersion
     else: dict['numpy'] = "None"

     # Check numpy
     (numpyVersion, numpyIncDir, numpyLibDir) = checkNumpy()
     if (numpyVersion != False): dict['numpy'] = numpyVersion
     else: dict['numpy'] = "None"
     
     # Check png
     (png, pngIncDir, pngLib) = checkPng(config.additionalLibPaths,
                                         config.additionalIncludePaths)
     if png: dict['png'] = pngLib
     else: dict['png'] = "None"

     # Check ffmpeg
     (mpeg, mpegIncDir, mpegLib) = checkMpeg(config.additionalLibPaths,
                                             config.additionalIncludePaths)
     if mpeg: dict['mpeg'] = mpegLib
     else: dict['mpeg'] = "None"
     
     # Check hdf5
     (hdf, hdfIncDir, hdfLib) = checkHdf(config.additionalLibPaths,
                                         config.additionalIncludePaths)
     if hdf: dict['hdf'] = hdfLib
     else: dict['hdf'] = "None"

     # Write dictionnary
     p.write("buildDict = "+str(dict))
     p.close()
     
#==============================================================================
# Ecrit la base d'installation (ancien config.py) dans le fichier
# installBase.py
# IN: dict: dictionnaire d'install
#==============================================================================
def writeInstallBase(dict):
     p = open("installBase.py", 'w')
     if p is None:
        raise SystemError("Error: can not open file installBase.py for writing.")

     # Write doc
     p.write("# This is the dictionary keeping track of installation.\n# The key is the machine name or ELSAPROD name. For each key a list is stored.\n# [description, \n# f77compiler, libfortdir, libfort, f90compiler, libf90dir, libf90, \n# Cppcompiler, libCpp, useOMP, \n# CPlotOffScreen \n# pngPath, mpegPath, adfPath, hdfPath].\n# Path are list of strings. useOMP, CPlotOffScreen are booleans. \n# Others are strings.\n")
    
     # Write dictionary
     #p.write("installDict = "+str(dict))

     # Pretty print dict
     p.write("installDict = {\n")
     keys = dict.keys()
     kc = 0
     for k in keys:
          p.write("###############################################################################\n")
          if isinstance(k, str): kstr = "\'%s\'"%k
          else: kstr = str(k)
          p.write("%s: [ "%kstr)
          list = dict[k]
          lc = 0
          for l in list:
               lc += 1
               if isinstance(l, str): lstr = "\'%s\'"%l
               else: lstr = str(l)
               if lc == 1:  p.write("%s,\n"%lstr)
               if lc == 2:  p.write("%s, # f77compiler\n"%lstr)
               elif lc == 3: p.write("%s, # f90compiler\n"%lstr)
               elif lc == 4: p.write("%s, # Cppcompiler\n"%lstr)
               elif lc == 5: p.write("%s, # CppAdditionalOptions\n"%lstr)
               elif lc == 6: p.write("%s, # f77AdditionalOptions\n"%lstr)
               elif lc == 7: p.write("%s, # useOMP\n"%lstr)
               elif lc == 8: p.write("%s, # static\n"%lstr)
               elif lc == 9: p.write("%s, # CPlotOffScreen\n"%lstr)
               elif lc == 10: p.write("%s, # additionalIncludePaths\n"%lstr)
               elif lc == 11: p.write("%s, # additionalLibs\n"%lstr)
               elif lc == 12: p.write("%s # additionalLibPaths\n"%lstr)
          kc += 1
          if kc == len(keys): p.write("]\n")
          else: p.write("], \n")
     p.write("}\n")
     p.close()

#==============================================================================
# Sur certains unix, le chemin d'installation contient un lib64
# Cree le symlink pour que lib et lib64 soit equivalent
#============================================================================== 
def symLinks():
     system = getSystem()[0]
     bits = getSystem()[1]
     if bits == '64': 
          try: import KCore.installPath as K
          except: import installPath as K
          libPath1 = K.libPath
          spls = libPath1.rsplit('/',1)
          if spls[1] == 'lib': libPath2 = spls[0]+'/lib64'
          elif spls[1] == 'lib64': libPath2 = spls[0]+'/lib'
          else: return
          ex1 = os.path.exists(libPath1)
          ex2 = os.path.exists(libPath2)
          lex1 = os.path.lexists(libPath1)
          lex2 = os.path.lexists(libPath2)
          if (not ex1 and lex1): # broken link 1
               os.remove(libPath1); ex1 = False
          if (not ex2 and lex2): # broken link 2
               os.remove(libPath2); ex2 = False
          if (ex1 and not ex2): os.symlink(libPath1, libPath2)
          elif (not ex1 and ex2): os.symlink(libPath2, libPath1)

#==============================================================================
# Cree les extensions d'un Module
# Cree une extension C/Api : Module.module (a partir de module.cpp)
# Des extensions Cython si necessaire: Module.f (a partir de la liste des pyx
# IN: module: module name ('Converter')
# IN: srcs: le module des sources
#==============================================================================
def createExtensions(module, srcs, includeDirs, libraryDirs, libraries, 
                     extraCompileArgs=[], extraLinkArgs=[]):
     from distutils.core import Extension
     listExtensions = []
     minor = module.lower()
     # C/Api module
     Extension(module+'.'+minor,
               sources=[module+'/'+minor+'.cpp'],
               include_dirs=[module]+includeDirs,
               library_dirs=libraryDirs,
               libraries=libraries,
               extra_compile_args=extraCompileArgs,
               extra_link_args=extraLinkArgs)
     # Cython extensions
     try: pyx_srcs = srcs.pyx_srcs
     except: pyx_srcs = []
     for i in pyx_srcs:
          f = i.replace('.pyx', '')
          Extension(module+'.'+f,
                    sources=[module+'/'+f+'.cpp'],
                    include_dirs=[module]+includeDirs,
                    library_dirs=libraryDirs,
                    libraries=libraries,
                    extra_compile_args=extraCompileArgs,
                    extra_link_args=extraLinkArgs)    
     return listExtensions

#==============================================================================
# Fortran builder
#==============================================================================
# Ajoute le fortran builder a env
def createFortranBuilder(env, dirs=[]):
     import SCons
     from SCons.Builder import Builder
     # Pre-processing
     path = ''
     for i in dirs: path += '"%s" -I'%i
     if path != '': path = path[:-3]
     bld = Builder(action=getPP()+'%s $SOURCES $TARGETS'%path, suffix='.f', 
                   src_suffix='.for')
     env.Append(BUILDERS={'FPROC': bld})
     # Fortran compiler
     fortran_builder = Builder(action='$FORTRANCOM',
                               suffix='.o', src_suffix='.f')
     env.Append(BUILDERS={'Fortran': fortran_builder})
     env.Replace(FORTRANFLAGS=getForArgs())
     env.Replace(FORTRANCOM='$FORTRAN $FORTRANFLAGS -c -o $TARGET $SOURCE')
     env.Replace(FORTRANSUFFIXES=['.f', '.F', '.f90', '.F90'])

     env.Replace(F90FLAGS=getForArgs())
     env.Replace(F95FLAGS=getForArgs())
     #env.Replace(SHF90=f90compiler)
     #env.Replace(SHF95=f90compiler)
     pref = getFortranModDirPrefix()
     env.Replace(FORTRANMODDIRPREFIX=pref) # doesnt work
     env.Replace(FORTRANMODDIR='MODS') # doesnt work
     args = getForArgs()
     if pref != '': args += [pref,'build'] # trick
     env.Replace(FORTRANFLAGS=args)
     return env

# Cree les noms des fichiers
def createFortranFiles(env, srcs):
     ppf = []
     try:
          for f in srcs.for_srcs:
               ffile = env.FPROC(target=f)
               ofile = env.Fortran(target=ffile)
               ppf.append(ofile[0])
     except: pass
     try:
          for f in srcs.f90_srcs:
               ofile = env.Fortran(target=f)
               ppf.append(ofile[0])
     except: pass
     return ppf

# Scan les .f pour faire les dependences (include)
def fortranScan(node, env, path, arg=None):
     import itertools
     import re
     # scan file to extract all possible include.
     contents = node.get_text_contents()
     names = [reo.findall(contents) for reo in [
               re.compile(r'^\s*from\s+(.+?)\s+cimport\s.*$', re.M),
               re.compile(r'^\s*cimport\s+(.+?)$', re.M),
               ]]
     names = itertools.chain(*names)
     # keep characters before " as ".
     names = [name.split(' as ')[0] for name in names]
     # split each cimport.
     names = itertools.chain(*[name.split(',') for name in names])
     names = [name.strip() for name in names]
     # remove duplications
     names = set(names)
     # prepend with the directory of the original pyx file.
     prefix = os.path.dirname(env.GetBuildPath(node))
     names = [os.path.join(prefix, '%s.pxd'%name) for name in names]
     # only take local pxd file and ignore anything unfound.
     names = [name for name in names if os.path.exists(name)]
     return [env.File(name) for name in names]

# Cree le scanner Fortran dans env
def createFortranScanner(env):
     import SCons
     fortranscanner = SCons.Scanner(function=fortranScan, skeys=['.f'], name='??')
     env.Append(SCANNERS=[fortranscanner])

#==============================================================================
# Builder Cython
#==============================================================================
def cythonSuffixEmitter(env, source):
     return "$CYTHONCFILESUFFIX"

# Cree le builder Cython dans env
def createCythonBuilder(env):
     import SCons
     from SCons.Builder import Builder
     from SCons.Action import Action
     incs = " -a --cplus "
     for i in env['CPPPATH']: incs += ' -I"%s" '%i

     cypath = ''
     if (env.has_key("CYTHONCOMPATH") and (env["CYTHONCOMPATH"] != "")):
          SYSPATH = ""
          if env.has_key("PYTHONPATH"): SYSPATH=env['PYTHONPATH']
          cypath = env["CYTHONCOMPATH"]
          if cypath != "": cypath += "/"
          pypath = "PYTHONPATH=%s:%s:%s "%(RESOURCELIBPATH,CYTHONLIBPATH,SYSPATH)
          cypath = pypath+cypath
     env["CYTHONCOM"] = cypath+"cython"+incs+" -o $TARGET $SOURCE"
     env["CYTHONCFILESUFFIX"] = ".cpp"

     c_file, cxx_file = SCons.Tool.createCFileBuilders(env)
     c_file.suffix['.pyx'] = cythonSuffixEmitter
     c_file.add_action('.pyx', Action("$CYTHONCOM"))

     cython = Builder(
          action=Action("$CYTHONCOM"),
          emitter={},
          suffix='.cpp', src_suffix='.pyx',
          single_source=1)
     env.Append(BUILDERS={'GenerateCython': cython})
     #env = createCythonScanner(env)
     return env

# Scan les .pyx pour faire les dependences (.pxd, .pxi)
def cythonScan(node, env, path, arg=None):
     import itertools
     import re
     # scan file to extract all possible cimports.
     contents = node.get_text_contents()
     names = [reo.findall(contents) for reo in [
               re.compile(r'^\s*from\s+(.+?)\s+cimport\s.*$', re.M),
               re.compile(r'^\s*cimport\s+(.+?)$', re.M),
               ]]
     names = itertools.chain(*names)
     # keep characters before " as ".
     names = [name.split(' as ')[0] for name in names]
     # split each cimport.
     names = itertools.chain(*[name.split(',') for name in names])
     names = [name.strip() for name in names]
     # remove duplications
     names = set(names)
     # prepend with the directory of the original pyx file.
     prefix = os.path.dirname(env.GetBuildPath(node))
     names = [os.path.join(prefix, '%s.pxd'%name) for name in names]
     # only take local pxd file and ignore anything unfound.
     names = [name for name in names if os.path.exists(name)]
     return [env.File(name) for name in names]

# Cree le scanner Cython dans env
def createCythonScanner(env):
     import SCons
     pyxscanner = SCons.Scanner(function=cythonScan, skeys=['.pyx'], name='PYX')
     env.Append(SCANNERS=[pyxscanner])

# Cree les targets des fichiers CPP
def createCythonFiles(env, srcs):
     cythonCpp = env.GenerateCython(srcs.pyx_srcs)
     deps = []
     for i in cythonCpp:
          base = os.path.dirname(str(i))
          deps += env.Install('../../'+base, i)
     return deps

#==============================================================================
if (__name__ == "__main__"):
   checkAll()
