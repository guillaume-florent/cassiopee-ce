# install les themes ttk
import KCore.installPath
import os
p = KCore.installPath.installPath+'/CPlot'
a = os.access(p+'/themes.tar', os.R_OK)
b = os.access(p+'/themes', os.R_OK)
if a and not b:
    os.system('cd \"%s\"; tar xvf themes.tar'%p)
