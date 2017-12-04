# - backgrounds -
import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}
VARS = []

#==============================================================================
# Si half = 1, cree que la moitie
# Si zpos = 1, position la boite en z min
# Sinon, cree tout
#==============================================================================
def createBox(half=0, zpos=0):
    alpha = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(alpha) == 1:
        ax = alpha[0]; ay = alpha[0]; az = alpha[0]
        bx = alpha[0]; by = alpha[0]; bz = alpha[0]
    elif len(alpha) == 3:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[0]; by = alpha[1]; bz = alpha[2]
    elif len(alpha) == 6:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[3]; by = alpha[4]; bz = alpha[5]
    else:
        CTK.TXT.insert('START', 'Borders factor incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    if zpos == 1: az = 0 # force zmin
    
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    LX = posEye[0]-posCam[0]
    LY = posEye[1]-posCam[1]
    LZ = posEye[2]-posCam[2]
    
    try: box = G.bbox(CTK.t)
    except: box = [0,0,0,1,1,1]
    hx = box[3]-box[0]
    hy = box[4]-box[1]
    hz = box[5]-box[2]
    if hx < 1.e-10: hx = 0.1
    if hy < 1.e-10: hy = 0.1
    if hz < 1.e-10: hz = 0.1
    hx = hx * 0.5; hy = hy * 0.5; hz = hz * 0.5
    ni = 3+ax+bx; nj = 3+ay+by; nk = 3+az+bz

    if LX < 0:
        # plan xmin
        b1 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (1,nj,nk) )
    else:
        # plan xmax
        b1 = G.cart( (box[3]+bx*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (1,nj,nk) )
    b1 = CPlot.addRender2Zone(b1, material='Solid',
                              color='White', meshOverlay=1)

    if LY < 0:
        # plan ymin
        b2 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,1,nk) )
    else:
        # plan ymax
        b2 = G.cart( (box[0]-ax*hx, box[4]+by*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,1,nk) )
    b2 = CPlot.addRender2Zone(b2, material='Solid',
                              color='White', meshOverlay=1)
    if LZ < 0:
        # plan zmin
        b3 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,nj,1) )
    else:
         b3 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[5]+bz*hz),
                      (hx,hy,hz), (ni,nj,1) )
    b3 = CPlot.addRender2Zone(b3, material='Solid',
                              color='White', meshOverlay=1)

    if half == 1: return [b1,b2,b3]

    if LX >= 0:
        # plan xmin
        b4 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (1,nj,nk) )
    else:
        # plan xmax
        b4 = G.cart( (box[3]+bx*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (1,nj,nk) )
    b4 = CPlot.addRender2Zone(b4, material='Solid',
                              color='White', meshOverlay=1)

    if LY >= 0:
        # plan ymin
        b5 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,1,nk) )
    else:
        # plan ymax
        b5 = G.cart( (box[0]-ax*hx, box[4]+by*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,1,nk) )
    b5 = CPlot.addRender2Zone(b5, material='Solid',
                              color='White', meshOverlay=1)
    if LZ >= 0:
        # plan zmin
        b6 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[2]-az*hz),
                     (hx,hy,hz), (ni,nj,1) )
    else:
        b6 = G.cart( (box[0]-ax*hx, box[1]-ay*hy, box[5]+bz*hz),
                     (hx,hy,hz), (ni,nj,1) )
    b6 = CPlot.addRender2Zone(b6, material='Solid',
                              color='White', meshOverlay=1)
    return [b1,b2,b3,b4,b5,b6]

#==============================================================================
def createZEllipse():
    import Geom.PyTree as D
    import Transform.PyTree as T
    alpha = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(alpha) == 1:
        ax = alpha[0]; ay = alpha[0]; az = alpha[0]
        bx = alpha[0]; by = alpha[0]; bz = alpha[0]
    elif len(alpha) == 3:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[0]; by = alpha[1]; bz = alpha[2]
    elif len(alpha) == 6:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[3]; by = alpha[4]; bz = alpha[5]
    else:
        CTK.TXT.insert('START', 'Borders factor incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    try: box = G.bbox(CTK.t)
    except: box = [0,0,0,1,1,1]
    hx = box[3]-box[0]
    hy = box[4]-box[1] 
    if hx < 1.e-10: hx = 0.1
    if hy < 1.e-10: hy = 0.1
  
    hx = hx * 0.7 * ax; hy = hy * 0.7 * ay
    cx = 0.5*(box[0]+box[3])
    cy = 0.5*(box[1]+box[4])
    if hx < hy:
        s = D.circle( (cx,cy,box[2]), hx, N=50)  
        s = T.contract(s, (cx,cy,box[2]), (1,0,0), (0,0,1), hy/hx ) 
    else:
        s = D.circle( (cx,cy,box[2]), hy, N=50)  
        s = T.contract(s, (cx,cy,box[2]), (0,1,0), (0,0,1), hx/hy )

    s = C.convertArray2Tetra(s); s = G.close(s, 1.e-6)
    p = G.fittingPlaster(s, bumpFactor=0.)
    s = G.gapfixer(s, p)
    s = CPlot.addRender2Zone(s, material='Solid',
                             color='White', meshOverlay=0)
    return [s]

#==============================================================================
def createZPlane():
    alpha = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(alpha) == 1:
        ax = alpha[0]; ay = alpha[0]; az = alpha[0]
        bx = alpha[0]; by = alpha[0]; bz = alpha[0]
    elif len(alpha) == 3:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[0]; by = alpha[1]; bz = alpha[2]
    elif len(alpha) == 6:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[3]; by = alpha[4]; bz = alpha[5]
    else:
        CTK.TXT.insert('START', 'Borders factor incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    try: box = G.bbox(CTK.t)
    except: box = [0,0,0,1,1,1]
    hx = box[3]-box[0]
    hy = box[4]-box[1] 
    if hx < 1.e-10: hx = 0.1
    if hy < 1.e-10: hy = 0.1
  
    hx = hx * 0.5; hy = hy * 0.5
    b = G.cart((box[0]-ax*hx, box[1]-ay*hy, box[2]), (hx,hy,1), (ax+bx+3,ay+by+3,1))
    b = CPlot.addRender2Zone(b, material='Solid',
                             color='White', meshOverlay=1)
    return [b]

#==============================================================================
def createGround():
    import Post.PyTree as P
    alpha = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(alpha) == 1:
        ax = alpha[0]; ay = alpha[0]; az = alpha[0]
        bx = alpha[0]; by = alpha[0]; bz = alpha[0]
    elif len(alpha) == 3:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[0]; by = alpha[1]; bz = alpha[2]
    elif len(alpha) == 6:
        ax = alpha[0]; ay = alpha[1]; az = alpha[2]
        bx = alpha[3]; by = alpha[4]; bz = alpha[5]
    else:
        CTK.TXT.insert('START', 'Borders factor incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    try: box = G.bbox(CTK.t)
    except: box = [0,0,0,1,1,1]
    hx = box[3]-box[0]
    hy = box[4]-box[1]
    hz = box[5]-box[2]
    if hx < 1.e-10: hx = 0.1
    if hy < 1.e-10: hy = 0.1
    if hz < 1.e-10: hz = 0.001
    h = max(hx, hy)
    ay = ax; bx = ax; by = ax; # force square
    deltax = 0.5*(h-hx); deltay = 0.5*(h-hy)

    hx = h * 0.5; hy = h * 0.5; hz = 0.1*hz
    b = G.cart((box[0]-ax*hx-deltax, box[1]-ay*hy-deltay, box[2]-hz), (hx,hy,hz), (ax+bx+3,ay+by+3,2))
    b = P.exteriorFaces(b)
    b = CPlot.addRender2Zone(b, material='Solid',
                             color='White', meshOverlay=1)
    # plafond (optionel)
    c = G.cart((box[0]-ax*hx-deltax, box[1]-ay*hy-deltay, box[2]+(ax+bx+2)*hx-hz), (hx,hy,hz), (ax+bx+3,ay+by+3,1))
    
    return [b, c]

#==============================================================================
def deleteBackgroundBase():
    nodes = Internal.getNodesFromName1(CTK.t, 'BACKGROUND')
    if nodes == []: return
    base = nodes[0]
    # delete from plotter
    zones = Internal.getNodesFromType(base, 'Zone_t')
    dels = []
    for z in zones: dels.append(base[0]+Internal.SEP1+z[0])
    CPlot.delete(dels)
    # delete from tree
    ret = Internal.getParentOfNode(CTK.t, base)
    del ret[0][2][ret[1]]
    
#==============================================================================
def setBackground(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    type = VARS[0].get()
    CTK.saveTree()
    if type == 'None':
        deleteBackgroundBase()
    else:
        deleteBackgroundBase()
        CTK.t = C.addBase2PyTree(CTK.t, 'BACKGROUND', 2)
        
        if type == 'Half-Box': B = createBox(1)
        elif type == 'Box': B = createBox(0)
        elif type == 'Z-Half-Box': B = createBox(1, 1)
        elif type == 'Z-Box': B = createBox(0, 1)
        elif type == 'Z-Ellipse': B = createZEllipse()
        elif type == 'Z-Plane': B = createZPlane()
        elif type == 'Z-Square-Ground': B = createGround()

        base = Internal.getNodesFromName1(CTK.t, 'BACKGROUND')[0]
        nob = C.getNobOfBase(base, CTK.t)
        for b in B: CTK.add(CTK.t, nob, -1, b)
        CTK.t = C.fillMissingVariables(CTK.t)
        
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkBackground', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Create a background.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Button-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkBackground')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Type de background -
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    if CTK.PREFS.has_key('tkBackgroundType'): 
        V.set(CTK.PREFS['tkBackgroundType'])
    # -1- Border
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    if CTK.PREFS.has_key('tkBackgroundBorder'): 
        V.set(CTK.PREFS['tkBackgroundBorder'])

    # - Type de background -
    B = TTK.OptionMenu(Frame, VARS[0], 'None', 'Half-Box', 'Box',
                       'Z-Half-Box', 'Z-Box', 'Z-Ellipse', 'Z-Plane',
                       'Z-Square-Ground',
                       command=setBackground)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of background.')

    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Borders N or (Ni;Nj;Nk) or (Ni-;Nj-;Nk-;Ni+;Nj+;Nk+).')
    B.bind('<Return>', setBackground)
    
    B = TTK.Button(Frame, text="Set", command=setBackground)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set or reset the background.')
    
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
    CTK.PREFS['tkBackgroundType'] = VARS[0].get()
    CTK.PREFS['tkBackgroundBorder'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('None')
    VARS[1].set('2')
    CTK.PREFS['tkBackgroundType'] = VARS[0].get()
    CTK.PREFS['tkBackgroundBorder'] = VARS[1].get()
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
    (win, menu, file, tools) = CTK.minimal('tkBackground '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
