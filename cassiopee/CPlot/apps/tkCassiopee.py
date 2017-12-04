# -- Cassiopee main app --
import Tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import os

# Liste des apps par sous menu et perso
TREEAPPS = ['tkTreeOps', 'tkCheckPyTree', '---',
            'tkFilter', '---',
            'tkFamily']
STATEAPPS = ['tkState', '---',
             'tkPrefs', 'tkPerfo', 'tkContainers', '---',
             #'tkLogFile', '---',
             'tkRuler', 'tkFind']
EDGEAPPS = ['tkCanvas', 'tkPoint', 'tkDraw','---',
            'tkExtractEdges', 'tkMapEdge']
SURFAPPS = ['tkBasicSurfs', 'tkText', '---',
            'tkFixer2', 'tkBoolean', 'tkSculpt', 'tkPaint', '---',
            'tkMapSurfs', 'tkFilterSurfs', 'tkSurfaceWalk', '---',
            'tkProjection']
MESHAPPS = ['tkCells', 'tkStretch', '---',
            'tkExtrusion', 'tkTetraMesher', 'tkTFI', 'tkSmooth', '---',
            'tkOctree', 'tkCollarMesh', 'tkBlader',
            #'tkPLM', 'tkPC1M',
            '---',
            'tkMeshQual', 'tkMeshInfo']
BLOCKAPPS = ['tkBlock', '---',
             'tkTransform', 'tkNGon', 'tkGhostCells', '---',
             'tkSplit', 'tkReorder']
BCAPPS = ['tkBC', '---',
          'tkChimera', 'tkIBC', '---',
          'tkExtractBC']
MOTIONAPPS = ['tkRigidMotion', 'tkTime']
SOLVERAPPS = ['tkInit', 'tkDistributor', 'tkDist2Walls', '---',
              'tkCassiopeeSolver', 'tkElsaSolver', 'tkFastSolver']
POSTAPPS = ['tkVariables', '---',
            'tkExtractMesh', '---',
            'tkStream', 'tkIsoLine', 'tkIsoSurf', '---',
            'tkInteg']
VISUAPPS = ['tkView', 'tkPlot', 'tkPlotXY', '---',
            'tkSlice', 'tkCellN', '---',
            'tkBackground']
RENDERAPPS = ['tkRenderSet', '---',
              'tkStereo', 'tkEffects', 'tkDemo', '---',
              'tkPovRay', 'tkLuxRender']

ALLAPPS = TREEAPPS + STATEAPPS + EDGEAPPS + SURFAPPS + MESHAPPS + \
          BLOCKAPPS + BCAPPS + MOTIONAPPS + SOLVERAPPS + POSTAPPS + \
          VISUAPPS + RENDERAPPS
PERSOAPPS = []

#==============================================================================
# Application global preference settings
#==============================================================================
def preferences():
    module = __import__('tkPrefs')
    module.createApp(win)
    module.activateApp()

#==============================================================================
# Add a personal app to pref file
#==============================================================================
def addPersonalApp():
    import tkFileDialog
    file = tkFileDialog.askopenfilename(
        filetypes=[('python', '*.py')])
    a = os.access(file, os.F_OK)
    if not a: return
    CTK.loadPrefFile()
    if CTK.PREFS.has_key('module'): CTK.PREFS['module'] += ' ;'+file
    else: CTK.PREFS['module'] = file
    CTK.savePrefFile()
    file = os.path.split(file)
    moduleName = os.path.splitext(file[1])[0]
    pathName = file[0]
    try:
        orig = sys.path; local = orig
        local.append(pathName)
        sys.path = local
        module = __import__(moduleName)
        CTK.TKMODULES[moduleName] = module
        sys.path = orig
        tools.add_command(label=moduleName,
                          command=module.showApp)
        module.createApp(frames[0])
        PERSOAPPS.append(moduleName)
        CTK.TXT.insert('START', 'Module %s added in tools menu\n'%moduleName)
    except:
        CTK.TXT.insert('START', 'Can not import '+moduleName+'.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
def notImplemented():
    CTK.TXT.insert('START', 'This functionality is not implemented.\n')
    CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# IN: app: nom d'applet (tkChose ou --- ou Sub/tkChose)
# si app=---, un separateur est affiche
# si app=Sub/tkChose, un sous-menu Sub est ajoute
# IN: menu: menu ou on ajoute l'app
# IN: frame: la frame de l'app
# IN: submenus: dict to keep trace of allready created submenus
# OUT: TKMODULES: le dictionnaire des modules importes
#==============================================================================
def addMenuItem(app, menu, frame, submenus):
    app = app.split('/')
    if len(app) == 2: submenu = app[0]; app = app[1]
    else: submenu = None; app = app[0]

    if submenu is None:
        if app == '---': menu.add_separator()
        else:
            #try:
            module = __import__(app)
            CTK.TKMODULES[app] = module
            name = app; name = '  '+name
            menu.add_command(label=name, command=module.showApp)
            module.createApp(frame)
            if auto[app] == 1: module.showApp()
            #except: pass
    else: # submenu
        if submenus.has_key(submenu):
            myMenu = submenus[submenu]
        else:
            myMenu = TK.Menu(menu, tearoff=0)
            submenus[submenu] = myMenu
            menu.add_cascade(label=submenu, menu=myMenu)
        if app == '---': myMenu.add_separator()
        else:
            #try:
            module = __import__(app)
            CTK.TKMODULES[app] = module
            name = app; name = '  '+name
            myMenu.add_command(label=name, command=module.showApp)
            module.createApp(frame)
            if auto[app] == 1: module.showApp()
            #except: pass

#==============================================================================
if (__name__ == "__main__"):

    # Ouverture du fichier de la ligne de commande
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE, density=1.)
            CTK.t = CTK.upgradeTree(CTK.t)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            fileName = os.path.split(CTK.FILE)[1]
            CPlot.CPlot.cplot.setWindowTitle(fileName)
        except: print 'Error: convertFile2PyTree: fail to read file %s.'%CTK.FILE

    # - Verifie l'arbre -
    errors = []
    if CTK.t != []:
        errors = Internal.checkPyTree(CTK.t, level=5)
        if errors == []: CTK.display(CTK.t)

    # Load and set prefs for interface
    CTK.loadPrefFile(); CTK.setPrefs()

    # Modules a ouvrir automatiquement
    auto = {}
    for app in ALLAPPS:
        app = app.split('/')
        if len(app) == 2: app = app[1]
        else: app = app[0]
        auto[app] = 0
    if CTK.PREFS.has_key('auto'):
        p = CTK.PREFS['auto']; p = p.split(';')
        for i in p:
            i = i.strip()
            auto[i] = 1

    # Main window
    (win, frames, menu, menus, file, tools) = CTK.minimal2('Cassiopee '+C.__version__,
                                                           show=False)

    # - Apps -
    submenus = {}
    for app in TREEAPPS: addMenuItem(app, menus[0], frames[0], submenus)
    submenus = {}
    for app in STATEAPPS: addMenuItem(app, menus[1], frames[1], submenus)
    submenus = {}
    for app in EDGEAPPS: addMenuItem(app, menus[2], frames[2], submenus)
    submenus = {}
    for app in SURFAPPS: addMenuItem(app, menus[3], frames[3], submenus)
    submenus = {}
    for app in MESHAPPS: addMenuItem(app, menus[4], frames[4], submenus)
    submenus = {}
    for app in BLOCKAPPS: addMenuItem(app, menus[5], frames[5], submenus)
    submenus = {}
    for app in BCAPPS: addMenuItem(app, menus[6], frames[6], submenus)
    submenus = {}
    for app in MOTIONAPPS: addMenuItem(app, menus[7], frames[7], submenus)
    submenus = {}
    for app in SOLVERAPPS: addMenuItem(app, menus[8], frames[8], submenus)
    submenus = {}
    for app in POSTAPPS: addMenuItem(app, menus[9], frames[9], submenus)
    submenus = {}
    for app in VISUAPPS: addMenuItem(app, menus[10], frames[10], submenus)
    submenus = {}
    for app in RENDERAPPS: addMenuItem(app, menus[11], frames[11], submenus)

    # Updated Apps from tree (containers from tree containers)
    if CTK.TKMODULES.has_key('tkContainers'): CTK.TKMODULES['tkContainers'].updateApp()

    # Get tkPlotXY if any
    CTK.TKPLOTXY = CTK.TKMODULES['tkPlotXY']

    # - Personal apps  -
    tools.add_command(label='Add a personal app',
                      command=addPersonalApp)

    import os.path, sys
    if CTK.PREFS.has_key('module'):
        mod = CTK.PREFS['module']
        mod = mod.split(';')
        for i in mod:
            i = i.strip()
            file = os.path.split(i)
            moduleName = os.path.splitext(file[1])[0]
            pathName = file[0]
            try:
                orig = sys.path; local = orig
                local.append(pathName)
                sys.path = local
                module = __import__(moduleName)
                CTK.TKMODULES[moduleName] = module
                sys.path = orig
                tools.add_command(label=moduleName,
                                  command=module.showApp)
                module.createApp(frames[0]); module.hideApp()
                PERSOAPPS.append(moduleName)
            except:
                CTK.TXT.insert('START', 'can not import '+moduleName+'.\n')
                CTK.TXT.insert('START', 'Error: ', 'Error')

    # Place win devant les autres fenetres
    win.deiconify(); win.focus_set()

    # - Erreur dans l'arbre -
    if errors != []:
        Panels.displayErrors(errors, header='Checking pyTree')
        CTK.t = Internal.correctPyTree(CTK.t, level=5)
        CTK.display(CTK.t)
        CTK.TKTREE.updateApp()

    # - Main loop -
    win.mainloop()

    # Del photos
    CTK.PHOTOS = []
