# - CPlot view settings -
import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Post.PyTree as P
import Transform.PyTree as T
import math

# local widgets list
WIDGETS = {}; VARS = []

# No de la variable courante
VARNO = -3
# Min/Max de cette variable pour tout l'arbre
VARMIN = 0.; VARMAX = 1.

#==============================================================================
def setMode(event=None):
    mode = VARS[6].get()
    imode = 0
    WIDGETS['mesh'].grid_forget()
    WIDGETS['solid'].grid_forget()
    WIDGETS['render'].grid_forget()
    WIDGETS['scalar'].grid_forget()
    WIDGETS['vector'].grid_forget()
    if mode == 'Mesh':
        imode = 0; WIDGETS['mesh'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'Solid':
        imode = 1; WIDGETS['solid'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'Render':
        imode = 2; WIDGETS['render'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'Scalar':
        imode = 3
        if VARS[18].get() == 'None':
            vars = C.getVarNames(CTK.t, excludeXYZ=True)
            if (len(vars)>0 and len(vars[0])>0): VARS[18].set(vars[0][0])
            displayField()
        WIDGETS['scalar'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'Vector':
        imode = 4; WIDGETS['vector'].grid(row=2, column=0, columnspan=3, sticky=TK.EW)
        if VARS[20].get() == 'None' or VARS[21].get() == 'None' or VARS[22].get() == 'None':
            vars = C.getVarNames(CTK.t, excludeXYZ=True)
            if len(vars)>0:
                lg = len(vars[0])
                if lg > 0: VARS[20].set(vars[0][0])
                if lg > 1: VARS[21].set(vars[0][1])
                elif lg > 0: VARS[21].set(vars[0][0])
                if lg > 2: VARS[22].set(vars[0][2])
                elif lg > 0: VARS[22].set(vars[0][0])
            displayVector()
    CPlot.setState(mode=imode)
    CTK.TXT.insert('START', 'Mode %s displayed.\n'%mode)

#==============================================================================
# Pour le champ scalaire (optionMenu)
def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    m = WIDGETS['scalarField'].children['menu']
    m.delete(0, TK.END)
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[18],
                      l=i:displayFieldl(v,l))

def updateVarNameList_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if (CTK.__MAINTREE__ <= 0 or nzs == []):
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    if WIDGETS.has_key('scalarField'):
        WIDGETS['scalarField']['values'] = allvars

#==============================================================================
def updateVarNameList__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if (CTK.__MAINTREE__ <= 0 or nzs == []):
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    m = WIDGETS['vectorField'+str(no)].children['menu']
    m.delete(0, TK.END)
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[19+no],
                      l=i:displayVectorl(v,l))

def updateVarNameList2__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if (CTK.__MAINTREE__ <= 0 or nzs == []):
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
   
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)

    if WIDGETS.has_key('vectorField'+str(no)):
        WIDGETS['vectorField'+str(no)]['values'] = allvars

#==============================================================================
def updateVarNameList1(event=None):
    updateVarNameList__(1)
def updateVarNameList2(event=None):
    updateVarNameList__(2)
def updateVarNameList3(event=None):
    updateVarNameList__(3)
def updateVarNameList1_2(event=None):
    updateVarNameList2__(1)
def updateVarNameList2_2(event=None):
    updateVarNameList2__(2)
def updateVarNameList3_2(event=None):
    updateVarNameList2__(3)

#==============================================================================
def displayFieldl(v, l):
    v.set(l); displayField()

#==============================================================================
def displayVectorl(v, l):
    v.set(l); displayVector()
    
#==============================================================================
# Display pour les scalar field
#==============================================================================
def displayField(event=None):
    if CTK.t == []: return
    global VARNO
    field = VARS[18].get()
    vars = C.getVarNames(CTK.t)[0]
    ifield = 0; lenvars = 0
    for i in vars:
        if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
            lenvars += 1
    for i in vars:
        if i == field: break
        if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
            ifield += 1
        
    if ifield == lenvars:
        CTK.TXT.insert('START', 'Variable not found in tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    VARNO = ifield
    compMin(); compMax()
    updateIsoWidgets()
    
    ifield = field.replace('centers:', '')
    CPlot.setState(mode=3, scalarField=ifield)
    CTK.TXT.insert('START', 'Variable %s displayed.\n'%field)

#==============================================================================
def displayVector(event=None):
    if CTK.t == []: return
    global VARNO
    field1 = VARS[20].get()
    field2 = VARS[21].get()
    field3 = VARS[22].get()
    vars = C.getVarNames(CTK.t)[0]
    ifield1 = 0; ifield2 = 0; ifield3 = 0
    lenvars = 0
    for i in vars:
        if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            lenvars += 1
            
    for i in vars:
        if i == field1: break
        if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            ifield1 += 1
        
    if (ifield1 == lenvars):
        CTK.TXT.insert('START', 'Variable not found in tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    for i in vars:
        if i == field2: break
        if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            ifield2 += 1
        
    if ifield2 == lenvars:
        CTK.TXT.insert('START', 'Variable not found in tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    for i in vars:
        if i == field3: break
        if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
            ifield3 += 1
        
    if ifield3 == lenvars:
        CTK.TXT.insert('START', 'Variable not found in tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    CPlot.setState(mode=4, vectorField1=ifield1, vectorField2=ifield2,
                   vectorField3=ifield3)
    CTK.TXT.insert('START', 'Variable '+field1+','+field2+','+field3+' displayed.\n')
    
#==============================================================================
def saveSlot():
    if CTK.t == []: return
    slot = VARS[0].get()
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    mode = VARS[6].get()
    field = VARS[18].get()
    niso = CTK.varsFromWidget(VARS[2].get(), type=2)
    niso = niso[0]
    isoEdges = WIDGETS['edges'].get()/25.-0.001
    light = VARS[5].get()
    if (light == 'IsoLight on'): light = 1
    else: light = 0
    colormap = VARS[4].get()
    CPlot._addRender2PyTree(CTK.t, slot=int(slot), posCam=posCam,
                            posEye=posEye, dirCam=dirCam,
                            mode=mode, scalarField=field, niso=niso,
                            isoEdges=isoEdges, isoLight=light,
                            colormap=colormap)
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Slot saved.\n')
     
#==============================================================================
def loadSlot():
    if CTK.t == []: return
    slot = VARS[0].get()
    sl = Internal.getNodesFromName2(CTK.t, 'Slot'+slot)
    if sl == []: return
    slot = sl[0]
    pos = Internal.getNodesFromName(slot, 'posCam')
    if pos != []:
        n = pos[0][1]; CPlot.setState(posCam=(n[0], n[1], n[2]))
    pos = Internal.getNodesFromName(slot, 'posEye')
    if pos != []:
        n = pos[0][1]; CPlot.setState(posEye=(n[0], n[1], n[2]))
    pos = Internal.getNodesFromName(slot, 'dirCam')
    if pos != []:
        n = pos[0][1]; CPlot.setState(dirCam=(n[0], n[1], n[2]))
    pos = Internal.getNodesFromName(slot, 'scalarField')
    if pos != []:
        n = pos[0][1]
        VARS[18].set(n.tostring())
    pos = Internal.getNodesFromName(slot, 'mode')
    if pos != []:
        n = pos[0][1]
        VARS[6].set(n.tostring())
        setMode()
    pos = Internal.getNodesFromName(slot, 'niso')
    if pos != []:
        n = pos[0][1]; niso = int(n[0])
        CPlot.setState(niso=niso)
        VARS[2].set(str(niso))
    pos = Internal.getNodesFromName(slot, 'isoEdges')
    if pos != []:
        n = pos[0][1]; CPlot.setState(isoEdges=n[0])
        WIDGETS['edges'].set( (n[0]+0.01)*50 )
    pos = Internal.getNodesFromName(slot, 'isoLight')
    if pos != []:
        n = pos[0][1]
        if (n[0] == 0): VARS[5].set('IsoLight off')
        else: VARS[5].set('IsoLight on')
    pos = Internal.getNodesFromName(slot, 'colormap')
    if pos != []:
        n = pos[0][1]
        VARS[4].set(n.tostring())
    setColormapLight()
    displayField()
    pos = Internal.getNodesFromName(slot, 'isoScales')
    if pos != []:
        updateIsoWidgets(); updateIsoPyTree()
    CTK.TXT.insert('START', 'Slot loaded.\n')

#==============================================================================
def setIsoLegend(event=None):
    if CTK.t == []: return
    legend = int(VARS[7].get())
    CPlot.setState(displayIsoLegend=legend)
   
#==============================================================================
def setNiso(event=None):
    if CTK.t == []: return
    niso = CTK.varsFromWidget(VARS[2].get(), type=2)
    if len(niso) != 1:
        CTK.TXT.insert('START', 'Number of solid isos is incorrect.\n'); return
    niso = niso[0]
    if VARNO >= 0: updateIsoPyTree()
    else: CPlot.setState(niso=niso)
    
#==============================================================================
def setIsoEdges(event=None):
    if CTK.t == []: return
    isoEdges = WIDGETS['edges'].get()/25.-0.001
    CPlot.setState(isoEdges=isoEdges)

#==============================================================================
def setColormapLight(event=None):
    colormap = VARS[4].get()
    light = VARS[5].get()
    style = 0
    if colormap == 'Blue2Red': style = 0
    elif colormap == 'Green2Red': style = 2
    elif colormap == 'Black2White': style = 4
    elif colormap == 'White2Black': style = 6
    elif colormap == 'Diverging': style = 8
    if light == 'IsoLight on': style += 1
    CPlot.setState(colormap=style)

#==============================================================================
def setScalarStyle(event=None):
    var = VARS[19].get()
    style = 0
    if var == 'Bands': style = 0
    elif var == 'Bands+mesh': style = 1
    elif var == 'Lines': style = 2
    elif var == 'Lines+mesh': style = 3
    CPlot.setState(scalarStyle=style)
    
#==============================================================================
def setDim(event=None):
    if CTK.t == []: return
    dim = CPlot.getState('dim')
    if (VARS[8].get() == '2D' and dim != 2): CPlot.setDim(2)
    if (VARS[8].get() == '3D' and dim != 3): CPlot.setDim(3)

#==============================================================================
def setXY(event=None):
    if CTK.t == []: return
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    dx = posCam[0]-posEye[0]
    dy = posCam[1]-posEye[1]
    dz = posCam[2]-posEye[2]
    d = math.sqrt(dx*dx + dy*dy + dz*dz)
    if dz > 0: d = -d
    posCam2 = (posEye[0], posEye[1] , posEye[2]+d)
    dirCam2 = (0,1,0)
    CPlot.setState(posCam=posCam2, dirCam=dirCam2)

#==============================================================================
def setYZ(event=None):
    if CTK.t == []: return
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    dx = posCam[0]-posEye[0]
    dy = posCam[1]-posEye[1]
    dz = posCam[2]-posEye[2]
    d = math.sqrt(dx*dx + dy*dy + dz*dz)
    if dx > 0: d = -d
    posCam2 = (posEye[0]+d, posEye[1] , posEye[2])
    dirCam2 = (0,0,1)
    CPlot.setState(posCam=posCam2, dirCam=dirCam2)

#==============================================================================
def setXZ(event=None):
    if CTK.t == []: return
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    dx = posCam[0]-posEye[0]
    dy = posCam[1]-posEye[1]
    dz = posCam[2]-posEye[2]
    d = math.sqrt(dx*dx + dy*dy + dz*dz)
    if dy < 0: d = -d
    posCam2 = (posEye[0], posEye[1]-d , posEye[2])
    dirCam2 = (0,0,1)
    CPlot.setState(posCam=posCam2, dirCam=dirCam2)

#==============================================================================
def compMin():
    global VARMIN
    if CTK.t == []: return
    if VARNO < 0: return
    var = VARS[18].get()
    VARMIN = C.getMinValue(CTK.t, var)
    
#==============================================================================
def compIsoMin(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    var = VARS[18].get()
    zones = CTK.getValidZones()
    if zones == []:
        if CTK.__MAINTREE__ == 1: 
            zones = Internal.getZones(CTK.t)
        else: zones = Internal.getZones(CTK.dt)
    if zones == []: return
    try:
        varmin = C.getMinValue(zones, var)
        VARS[9].set(str(varmin))
        delta = max(VARMAX-VARMIN, 1.e-12)
        s = 100*(varmin-VARMIN)/delta
        WIDGETS['min'].set(s)
        updateIsoPyTree()
    except: pass

#==============================================================================
def compMax():
    global VARMAX
    if CTK.t == []: return
    if VARNO < 0: return
    var = VARS[18].get()
    VARMAX = C.getMaxValue(CTK.t, var)
    
#==============================================================================
def compIsoMax(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    var = VARS[18].get()
    zones = CTK.getValidZones()
    if zones == []:
        if CTK.__MAINTREE__ == 1: 
            zones = Internal.getZones(CTK.t)
        else: zones = Internal.getZones(CTK.dt)
    if zones == []: return
    try:
        varmax = C.getMaxValue(zones, var)
        VARS[10].set(str(varmax))
        delta = max(VARMAX-VARMIN, 1.e-12)
        s = 100*(varmax-VARMIN)/delta
        WIDGETS['max'].set(s)
        updateIsoPyTree()
    except: pass
    
#==============================================================================
def setIsoMin(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    updateIsoPyTree()
    
#==============================================================================
def setIsoMax(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    updateIsoPyTree()

#==============================================================================
def scaleIsoMin(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    val = WIDGETS['min'].get()
    fmin = VARMIN + (VARMAX - VARMIN)*val/100.
    VARS[9].set(str(fmin))
    updateIsoPyTree()
    
#==============================================================================
def scaleIsoMax(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    val = WIDGETS['max'].get()
    fmax = VARMIN + (VARMAX - VARMIN)*val/100.
    VARS[10].set(str(fmax))
    updateIsoPyTree()
    
#==============================================================================
# update iso widgets from pyTree and VARNO
def updateIsoWidgets():
    slot = VARS[0].get()
    sl = Internal.getNodeFromName2(CTK.t, 'Slot'+slot)
    if sl is None:
        compIsoMin(); compIsoMax(); return
    pos = Internal.getNodeFromName(sl, 'isoScales')
    if pos is not None:
        n = pos[1]; l = n.shape[0]
        list = []; c = 0
        while c < l:
            if n[c] == VARNO:
                VARS[2].set(str(int(n[c+1])))
                VARS[9].set(str(n[c+2]))
                VARS[10].set(str(n[c+3]))
                delta = max(VARMAX-VARMIN, 1.e-12)
                s = 100*(n[c+2]-VARMIN)/delta
                WIDGETS['min'].set(s)
                s = 100*(n[c+3]-VARMIN)/delta
                WIDGETS['max'].set(s)
                return
            c += 4
        compIsoMin(); compIsoMax()
    else:
        compIsoMin(); compIsoMax()

#==============================================================================
# update pyTree from iso widgets; display; update tree
def updateIsoPyTree():
    if VARNO >= 0:
        niso = int(VARS[2].get())
        fmin = float(VARS[9].get())
        fmax = float(VARS[10].get())
        list = [VARNO, niso, fmin, fmax]
        CPlot.setState(isoScales=list)
        slot = int(VARS[0].get())
        CPlot._addRender2PyTree(CTK.t, slot=slot, isoScales=list)
        CTK.TKTREE.updateApp()

#==============================================================================
def deleteIsoBase():
    nodes = Internal.getNodesFromName1(CTK.t, 'ISOLINES')
    if nodes == []: return
    base = nodes[0]
    # delete from plotter
    zones = Internal.getNodesFromType1(base, 'Zone_t')
    dels = []
    for z in zones: dels.append('%s/%s'%(base[0],z[0]))
    CPlot.delete(dels)
    ret = Internal.getParentOfNode(CTK.t, base)
    del ret[0][2][ret[1]]
    
#==============================================================================
def triggerIsoLines(event=None):
    if CTK.t == []: return
    if VARNO < 0: return
    isoLines = VARS[11].get()
    if isoLines == '0':
        deleteIsoBase()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    else:
        deleteIsoBase()
        field = VARS[18].get()
        nlevels = int(VARS[2].get())+1
        fmin = float(VARS[9].get())
        fmax = float(VARS[10].get())
        zones = CTK.getValidZones()
        isos = []
        for v in xrange(nlevels):
            value = fmin + (fmax-fmin)/(nlevels-1)*v
            for zone in zones:
                try:
                    i = P.isoLine(zone, field, value)
                    isos.append(i)
                except: pass
        CTK.t = C.addBase2PyTree(CTK.t, 'ISOLINES', 1)
        bases = Internal.getNodesFromName1(CTK.t, 'ISOLINES')
        if isos != []:
            isos = T.join(isos)
            base = bases[0]
            nob = C.getNobOfBase(base, CTK.t)
            CTK.add(CTK.t, nob, -1, isos)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CPlot.render()

#==============================================================================
def setVals():
    if VARS[12].get() == '1': CPlot.setState(edgifyActivatedZones=1)
    else: CPlot.setState(edgifyActivatedZones=0)
    if VARS[13].get() == '1': CPlot.setState(edgifyDeactivatedZones=1)
    else: CPlot.setState(edgifyDeactivatedZones=0)

#==============================================================================
def setStyle(event=None):
    style = 0; v = VARS[16].get()
    if v == 'Monocolor wires+solid': style = 0
    elif v == 'Multicolor wireframes': style = 1
    elif v == 'Multicolor wires+solid': style = 2
    elif v == 'Black wires+solid': style = 3
    elif v == 'White wires+solid': style = 4
    CPlot.setState(meshStyle=style)
    style = 0; v = VARS[17].get()
    if v == 'Monocolor/1-side': style = 0
    elif v == 'Multicolor/2-sides': style = 1
    elif v == 'White/2-sides': style = 3
    elif v == 'Black/2-sides': style = 4
    CPlot.setState(solidStyle=style)
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()
    
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkView', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage view.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Button-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkView')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Slot -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -1- Axis -
    V = TK.StringVar(win); V.set('Default'); VARS.append(V)
    # -2- Niso
    V = TK.StringVar(win); V.set('25'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewNiso'): V.set(CTK.PREFS['tkViewNiso'])
    # -3- isoEdges
    V = TK.StringVar(win); V.set('-0.5'); VARS.append(V)
    # -4- colormap type
    V = TK.StringVar(win); V.set('Blue2Red'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewColormap'): V.set(CTK.PREFS['tkViewColormap'])
    # -5- iso light
    V = TK.StringVar(win); V.set('IsoLight on'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewIsoLight'): V.set(CTK.PREFS['tkViewIsoLight'])
    # -6- Displayed mode
    V = TK.StringVar(win); V.set('Mesh'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewMode'): V.set(CTK.PREFS['tkViewMode'])
    # -7- Legende
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewLegend'): V.set(CTK.PREFS['tkViewLegend'])
    # -8- Dim
    V = TK.StringVar(win); V.set('3D'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewDim'): V.set(CTK.PREFS['tkViewDim'])
    # -9- Min value of isos
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewMin'): V.set(CTK.PREFS['tkViewMin'])
    # -10- Max value of isos
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewMax'): V.set(CTK.PREFS['tkViewMax'])
    # -11- IsoLine trigger
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -12- Edge for activated zones
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewEdgeA'): V.set(CTK.PREFS['tkViewEdgeA'])
    # -13- Edge for deactivated zones
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewEdgeD'): V.set(CTK.PREFS['tkViewEdgeD'])
    # -14- Shadow
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -15- DOF
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -16- Mesh style
    V = TK.StringVar(win); V.set('Multicolor wires+solid'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewMeshStyle'): 
        V.set(CTK.PREFS['tkViewMeshStyle'])
    # -17- Solid style
    V = TK.StringVar(win); V.set('Monocolor/1-side'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewSolidStyle'): 
        V.set(CTK.PREFS['tkViewSolidStyle'])
    # -18- scalar variable
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    # -19- Scalar style
    V = TK.StringVar(win); V.set('Bands'); VARS.append(V)
    if CTK.PREFS.has_key('tkViewScalarStyle'): 
        V.set(CTK.PREFS['tkViewScalarStyle'])
    # -20- vector variable 1
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    # -21- vector variable 2
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    # -22- vector variable 3
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    
    # - Dimension '2D ou 3D'
    B = TTK.OptionMenu(Frame, VARS[8], '3D', '2D',
                       command=setDim)
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)

    # - Choix du mode -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    B = TTK.OptionMenu(F, VARS[6], 'Mesh', 'Solid', 'Render',
                       'Scalar', 'Vector', command=setMode)
    B.grid(sticky=TK.EW)
    F.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Displayed mode.')
    WIDGETS['mode'] = B

    # - Axis -
    B = TTK.Button(Frame, text="XY", command=setXY)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(Frame, text="YZ", command=setYZ)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(Frame, text="XZ", command=setXZ)
    B.grid(row=0, column=2, sticky=TK.EW)

    # - Mesh frame -
    Mesh = TTK.LabelFrame(Frame, borderwidth=0)
    Mesh.columnconfigure(0, weight=0)
    Mesh.columnconfigure(1, weight=1)
    Mesh.grid(row=2, column=0, columnspan=3)
    WIDGETS['mesh'] = Mesh

    # - Mesh settings -
    B = TTK.Label(Mesh, text="Style:")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Mesh style.')
    B = TTK.OptionMenu(Mesh, VARS[16], 'Monocolor wires+solid', 'Multicolor wireframes', 'Multicolor wires+solid', 'Black wires+solid', 'White wires+solid', command=setStyle)
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Solid frame -
    Solid = TTK.LabelFrame(Frame, borderwidth=0)
    Solid.columnconfigure(0, weight=0)
    Solid.columnconfigure(1, weight=1)
    WIDGETS['solid'] = Solid

    # - Solid settings -
    B = TTK.Label(Solid, text="Style:")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Solid style.')
    B = TTK.OptionMenu(Solid, VARS[17], 'Monocolor/1-side', 'Multicolor/2-sides', 'White/2-sides', 'Black/2-sides', command=setStyle)
    B.grid(row=0, column=1, sticky=TK.EW)
    
    # - Render frame -
    Render = TTK.LabelFrame(Frame, borderwidth=0)
    Render.columnconfigure(0, weight=1)
    Render.columnconfigure(1, weight=1)
    WIDGETS['render'] = Render

    # - Render settings -
    B = TTK.Button(Render, text="Render panel", command=Panels.openRenderPanel)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Open render panel.')

    # - Scalar field frame -
    Scalar = TTK.LabelFrame(Frame, borderwidth=0)
    Scalar.columnconfigure(0, weight=0)
    Scalar.columnconfigure(1, weight=1)
    WIDGETS['scalar'] = Scalar

    # - Scalar field settings -
    B = TTK.Label(Scalar, text='Field:')
    B.grid(row=0, column=0, sticky=TK.EW)
    F = TTK.Frame(Scalar, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[18], '', command=displayField)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['scalarField'] = B
    else:
       B = ttk.Combobox(F, textvariable=VARS[18], 
                        values=[], state='readonly', width=10)
       B.bind('<<ComboboxSelected>>', displayField) 
       B.grid(sticky=TK.EW)
       F.bind('<Enter>', updateVarNameList_2)
       F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
       WIDGETS['scalarField'] = B

    B = TTK.Label(Scalar, text='Style:')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Scalar, VARS[19], 'Bands', 'Bands+mesh', 
                       'Lines', 'Lines+mesh',
                       command=setScalarStyle)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Scalar style.')
    
    B = TTK.Label(Scalar, text='Isos:')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Scalar, textvariable=VARS[2], width=4, background='White')
    B.bind('<Return>', setNiso)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of solid isos.')
    #B = TK.Entry(Scalar, textvariable=VARS[3], width=4, background='White')
    #B.bind('<Return>', setIsoEdges)
    #B.grid(row=2, column=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Width of edges (-0.5:deactivate, 1.:ok).')
    B = TTK.Scale(Scalar, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setIsoEdges, value=0)
    WIDGETS['edges'] = B
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Width of isolines.')
    
    B = TTK.Button(Scalar, text="Min", command=compIsoMin)
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Entry(Scalar, textvariable=VARS[9], width=4, background='White')
    B.bind('<Return>', setIsoMin)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Min of isos for this field.')
    B = TTK.Scale(Scalar, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=scaleIsoMin, showvalue=0, borderwidth=1, value=0)
    WIDGETS['min'] = B
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Min of isos for this field.')
    
    B = TTK.Button(Scalar, text="Max", command=compIsoMax)
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(Scalar, textvariable=VARS[10], width=4, background='White')
    B.bind('<Return>', setIsoMax)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max of isos for this field.')
    B = TTK.Scale(Scalar, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=scaleIsoMax, showvalue=0, value=100)
    WIDGETS['max'] = B
    B.grid(row=4, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max of isos for this field.')
    
    # - Colormap + light -
    B = TTK.Checkbutton(Scalar, text='Leg', variable=VARS[7],
                        command=setIsoLegend)
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Display color legend.')
    B = TTK.OptionMenu(Scalar, VARS[4], 'Blue2Red', 'Green2Red',
                       'Black2White', 'White2Black', 'Diverging',
                       command=setColormapLight)
    B.grid(row=5, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Colormap type.')
    B = TTK.OptionMenu(Scalar, VARS[5], 'IsoLight off', 'IsoLight on',
                       command=setColormapLight)
    B.grid(row=5, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='IsoLight.')

    # - Vector field frame -
    Vector = TTK.LabelFrame(Frame, borderwidth=0)
    Vector.columnconfigure(0, weight=0)
    Vector.columnconfigure(1, weight=1)
    WIDGETS['vector'] = Vector

    # - Vector field settings -
    B = TTK.Label(Vector, text='Field1:')
    B.grid(row=0, column=0, sticky=TK.EW)
    F = TTK.Frame(Vector, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[20], '', command=displayVector)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField1'] = B
        B = TK.Label(Vector, text='Field2:')
        B.grid(row=1, column=0, sticky=TK.EW)
        F = TK.Frame(Vector, borderwidth=0)
        F.columnconfigure(0, weight=1)
        B = TK.OptionMenu(F, VARS[21], '', command=displayVector)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField2'] = B
        B = TK.Label(Vector, text='Field3:')
        B.grid(row=2, column=0, sticky=TK.EW)
        F = TK.Frame(Vector, borderwidth=0)
        F.columnconfigure(0, weight=1)
        B = TK.OptionMenu(F, VARS[22], '', command=displayVector)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3)
        F.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField3'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[20], 
                         values=[], state='readonly', width=15)
        B.bind('<<ComboboxSelected>>', displayVector) 
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1_2)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField1'] = B
        B = TK.Label(Vector, text='Field2:')
        B.grid(row=1, column=0, sticky=TK.EW)
        F = TK.Frame(Vector, borderwidth=0)
        F.columnconfigure(0, weight=1)
        B = ttk.Combobox(F, textvariable=VARS[21], 
                         values=[], state='readonly', width=15)
        B.bind('<<ComboboxSelected>>', displayVector) 
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2_2)
        F.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField2'] = B
        B = TK.Label(Vector, text='Field3:')
        B.grid(row=2, column=0, sticky=TK.EW)
        F = TK.Frame(Vector, borderwidth=0)
        F.columnconfigure(0, weight=1)
        B = ttk.Combobox(F, textvariable=VARS[22], 
                         values=[], state='readonly', width=15)
        B.bind('<<ComboboxSelected>>', displayVector) 
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3_2)
        F.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['vectorField3'] = B

    # - Edge activation -
    B = TTK.Checkbutton(Frame, text='AEdges', variable=VARS[12],
                        command=setVals)
    BB = CTK.infoBulle(parent=B, text='Show edges for activated zones.')
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Checkbutton(Frame, text='DEdges', variable=VARS[13],
                        command=setVals)
    BB = CTK.infoBulle(parent=B, text='Show edges for deactivated zones.')
    B.grid(row=3, column=1, sticky=TK.EW)

    # - Slot -
    B = TTK.OptionMenu(Frame, VARS[0], '0', '1', '2', '3', '4', '5')
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Slot.')
    B = TTK.Button(Frame, text="Save", command=saveSlot)
    B.grid(row=5, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Save view, ... in slot.')
    B = TTK.Button(Frame, text="Load", command=loadSlot)
    B.grid(row=5, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Load view, ... from slot.')

    if CTK.PREFS.has_key('tkViewMode'): setMode()

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkViewNiso'] = VARS[2].get()
    CTK.PREFS['tkViewColormap'] = VARS[4].get()
    CTK.PREFS['tkViewIsoLight'] = VARS[5].get()
    CTK.PREFS['tkViewMode'] = VARS[6].get()
    CTK.PREFS['tkViewLegend'] = VARS[7].get()
    CTK.PREFS['tkViewDim'] = VARS[8].get()
    CTK.PREFS['tkViewMin'] = VARS[9].get()
    CTK.PREFS['tkViewMax'] = VARS[10].get()
    CTK.PREFS['tkViewEdgeA'] = VARS[12].get()
    CTK.PREFS['tkViewEdgeD'] = VARS[13].get()
    CTK.PREFS['tkViewMeshStyle'] = VARS[16].get()
    CTK.PREFS['tkViewSolidStyle'] = VARS[17].get()
    CTK.PREFS['tkViewScalarStyle'] = VARS[19].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[2].set('25')
    VARS[4].set('Blue2Red')
    VARS[5].set('IsoLight on')
    VARS[6].set('Mesh')
    VARS[7].set('0')
    VARS[8].set('3D')
    VARS[9].set('0.')
    VARS[10].set('1.')
    VARS[12].set('0')
    VARS[13].set('0')
    VARS[16].set('Monocolor wires+solid')
    VARS[17].set('Monocolor/1-side')
    VARS[19].set('Bands')
    CTK.PREFS['tkViewNiso'] = VARS[2].get()
    CTK.PREFS['tkViewColormap'] = VARS[4].get()
    CTK.PREFS['tkViewIsoLight'] = VARS[5].get()
    CTK.PREFS['tkViewMode'] = VARS[6].get()
    CTK.PREFS['tkViewLegend'] = VARS[7].get()
    CTK.PREFS['tkViewDim'] = VARS[8].get()
    CTK.PREFS['tkViewMin'] = VARS[9].get()
    CTK.PREFS['tkViewMax'] = VARS[10].get()
    CTK.PREFS['tkViewEdgeA'] = VARS[12].get()
    CTK.PREFS['tkViewEdgeD'] = VARS[13].get()
    CTK.PREFS['tkViewMeshStyle'] = VARS[16].get()
    CTK.PREFS['tkViewSolidStyle'] = VARS[17].get()
    CTK.PREFS['tkViewScalarStyle'] = VARS[19].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if (__name__ == "__main__"):
    import sys
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkView '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
