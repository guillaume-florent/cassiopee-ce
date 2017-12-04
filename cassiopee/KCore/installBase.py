# This is the dictionary keeping track of installation.
# The key is the machine name. For each key a list is stored.
# [description,
# f77compiler, f90compiler, Cppcompiler, useOMP, static, CPlotOffScreen,
# additionalIncludePaths, additionalLibs, additionalLibPaths].
# Paths are list of strings. useOMP, CPlotOffScreen are booleans.
# Others are strings.
installDict = {
###############################################################################
'WDSNA81OZ': [ 'Machine de production win32',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['c:/MinGW/include'], # additionalIncludePaths
['gfortran', 'gomp', 'pthread'], # additionalLibs
['c:/MinGW/lib', 'c:/Python27/libs'] # additionalLibPaths
],
###############################################################################
'node6.cluster': [ 'MacOSX',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
],
###############################################################################
'd1cn0001':[ 'AIRBUS HPC4B dev/val cluster',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
False, # useOMP
False, # static
False, # CPlotOffScreen
['/opt/hpmpi/include', '/opt/soft/cdtng/tools/portage/1.10.1_3D/usr/include'], # additionalIncludePaths
['svml', 'irc', 'ifcore'], # additionalLibs
['/opt/soft/cdtng/tools/intelcompiler/14.0/compiler/lib/intel64', '/opt/soft/cdtng/tools/portage/1.10.1_3D/usr/lib', '/opt/hpmpi/lib/linux_amd64'] # additionalLibPaths
],
###############################################################################
'caefr0p045': [ 'AIRBUS GISEH',
'ifort', # f77compiler
'ifort', # libf90dir
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/opt/soft/cdtng/tools/portage/1.10.1_3D/usr/include'], # additionalIncludePaths
['irc', 'mpi'], # additionalLibs
['/softs/compilers/intel/11.0/084/lib/intel64','/opt/soft/cdtng/tools/portage/1.10.1_3D/usr/lib'] # additionalLibPaths
],
###############################################################################
'wfrontend1': [ 'Safran Kairos',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
False, # useOMP
False, # static
False, # CPlotOffScreen
['/appl/APPLI_SNECMA/HDF5/oper/1.8.11/include'], # additionalIncludePaths
['ifcore', 'svml', 'irc'], # additionalLibs
['/opt/intel/composer_xe_2013_sp1.0.080/lib/intel64', '/appl/APPLI_SNECMA/HDF5/oper/1.8.11/lib'] # additionalLibPaths
],
###############################################################################
'santafe': [ 'MacOSX',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/X11/include' ], # additionalIncludePaths
['python2.7', 'ifcore'], # additionalLibs
['/usr/X11/lib', '/System/Library/Frameworks/OpenGL.framework/Libraries/'], # additionalLibPaths
],
###############################################################################
'daapuv': [ 'Onera DAAP',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/usr/local/hdf5-1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.7/lib'] # additionalLibPaths
],
###############################################################################
'celeste': [ 'Grosse machine de post-traitement onera',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPath
],
###############################################################################
'oneroa142': [ 'Onera',
'/opt/intel/fc/9.1.036/bin/ifort', # f77compiler
'/opt/intel/fc/9.1.036/bin/ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'linux64': [ 'Production linux64',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/stck1/benoit/include'], # additionalIncludePaths
[], # additionalLibs
['/stck1/benoit/lib'] # additionalLibPaths
],
###############################################################################
'eos...z': [ 'Onera-eosXXXz',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/local/hdf5-gnu-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-gnu-1.8.8/lib'] # additionalLibPaths
],
###############################################################################
'eos...': [ 'Onera-eos (legacy-doit etre apres eosZ)',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'ld...': [ 'Onera-ld',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'tiamat': [ 'Onera machine elsA',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=16'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/home/benoit/x86_64t/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/x86_64t'] # additionalLibPaths
],
###############################################################################
'austri.onera': [ 'Onera machine austri',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/home/benoit/aus/include'], # additionalIncludePaths
[], # additionalLibs
['/home/benoit/aus/lib'] # additionalLibPaths
],
###############################################################################
'westri': [ 'Onera machine westri',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=64'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/usr/local/hdf5-1.8.8-intel-16/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.8-intel-16/lib'] # additionalLibPaths
],
###############################################################################
'giulia': [ 'Onera machine elsA-ASO',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
["/tmp_opt/lib/hdf5-1.8.8-intel-16-impi/include",
"/usr/local/intel/studio/2016/compilers_and_libraries_2016.0.109/linux/mpi/include64"], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'mangrove': [ 'Onera machine GPU',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'moloch': [ 'Onera machine Cedre',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'cc-wdsna': [ 'Portable redhat',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[], # additionalLibPaths
],
###############################################################################
'cephee': [ 'Onera Cassiopee cluster',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/home/tools/local/x86_64a/include'], # additionalIncludePaths
[], # additionalLibs
['/home/tools/local/x86_64a/lib'] # additionalLibPaths
],
###############################################################################
'btmclx2': [ 'Turbomeca cluster',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
['/usr/lib/gcc/x86_64-redhat-linux/4.1.2'] # additionalLibPaths
],
###############################################################################
'WDSNAXXX': [ '??',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'visio': [ 'Onera gfx',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
['Xxf86vm'], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'WDSNA917Z': [ 'Machine de production win64',
'x86_64-w64-mingw32-gfortran', # f77compiler
'x86_64-w64-mingw32-gfortran', # f90compiler
'x86_64-w64-mingw32-gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
True, # static
False, # CPlotOffScreen
['c:/TDM-GCC-64/include'], # additionalIncludePaths
['gfortran', 'gomp', 'quadmath'], # additionalLibs
['c:/TDM-GCC-64/lib', 'c:/Python2.7/libs'] # additionalLibPaths
#['c:/TDM-GCC-64/lib', 'c:/Users/Adminstrateur/Anaconda2/libs'] # additionalLibPaths
],
###############################################################################
'fulvio': [ 'Onera gfx',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/local/hdf5-intel-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-intel-1.8.8/lib'] # additionalLibPaths
],
###############################################################################
'ubuntu': [ 'Generic ubuntu',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'linux': [ 'Generic linux',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'pdev': [ 'Machine Airbus',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/opt/soft/cdtng/tools/portage/1.9/usr/include', '/opt/hpmpi/include'], # additionalIncludePaths
[], # additionalLibs
['/opt/soft/cdtng/tools/portage/1.9/usr/lib', '/opt/hpmpi/lib', '/opt/soft/cdtng/tools/intelcompiler/11.0/lib/intel64'] # additionalLibPaths
],
###############################################################################
'laura': [ 'Machine acou',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/local/hdf5/1.8.7/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/lib64', '/usr/local/hdf5/1.8.7/lib','/tmp_opt/Python/2.7.3/icc-mpt/lib'] # additionalLibPaths
],
###############################################################################
'service': [ 'Onera Stelvio',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
["-axAVX,SSE4.2"], # CppAdditionalOptions
["-axAVX,SSE4.2"], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
'r.i.n.': [ 'Onera Stelvio-batch node',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
["-axAVX,SSE4.2"], # CppAdditionalOptions
["-axAVX,SSE4.2"], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'sator': [ 'Onera machine Sator',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
['-DCACHELINE=32','-Dvtune'], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
['/tmp_opt/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/include/','/tmp_opt/intel/studio/2017/advisor_2017.1.3.510716/include/intel64'], # additionalIncludePaths
['ittnotify','advisor'], # additionalLibs
['/tmp_opt/intel/studio/2017/vtune_amplifier_xe_2017.3.0.510739/lib64/','/tmp_opt/intel/studio/2017/advisor_2017.1.3.510716/lib64','/tmp_opt/lib/hdf5-1.8.17-intel-17/lib/'] # additionalLibPaths
], 
###############################################################################
'chi85bi': [ 'Machine EDF',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'Raspail': [ 'Onera DTIM',
'gfortran-7', # f77compiler
'gfortran-7', # f90compiler
#'clang++-4.0', # Cppcompiler
'g++-7', # Cppcompiler
['-std=c++98', '-pedantic', '-march=native'], # CppAdditionalOptions
['-march=native', '-fdefault-real-8', '-fdefault-double-8'], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/include/hdf5/serial/'], # additionalIncludePaths
['gfortran','gomp'], # additionalLibs
['/usr/lib/gcc/x86_64-linux-gnu/7',
 '/usr/lib/x86_64-linux-gnu/'] # additionalLibPaths
],
###############################################################################
'curie': [ 'Curie',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/usr/local/hdf5-1.8.8/include'], # additionalIncludePaths
[], # additionalLibs
['/usr/local/hdf5-1.8.8/lib'] # hdfPath
],
 ##############################################################################
'madmax64': [ 'Onera cluster madmax',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
True, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
['/usr/local/intel/cluster_studio/2012_0_032/lib/intel64'], # additionalLibPaths
],
###############################################################################
'localhost.localdomain': [ 'Unknown',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'default': [ 'Default',
'gfortran', # f77compiler
'gfortran', # f90compiler
'gcc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
[], # additionalIncludePaths
[], # additionalLibs
[] # additionalLibPaths
],
###############################################################################
'stelvio_impi15': [ 'Onera Stelvio Full intel',
'ifort', # f77compiler
'ifort', # f90compiler
'icc', # Cppcompiler
[], # CppAdditionalOptions
[], # f77AdditionalOptions
True, # useOMP
False, # static
False, # CPlotOffScreen
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/include',
 # '/tmp_opt/lib/hdf5/1.8.17/15/impi/include',
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/include'], # additionalIncludePaths
['mpi'], # additionalLibs
['/tmp_opt/lib/hdf5-1.8.8-intel-15-impi/lib', 
 '/tmp_opt/lib/hdf5/1.8.17/15/impi/lib', 
 '/tmp_opt/intel/studio/2015/impi/5.0.3.048/intel64/lib'] # additionalLibPaths
]
}
