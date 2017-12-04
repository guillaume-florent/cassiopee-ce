# - check pyTree integrity -
import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}
VARS = []
    
#==============================================================================
def runCheckPyTree():
    if CTK.t == []: return
    errors = []
    v = VARS[3].get()
    if (v == 'All conformity' or v == ' > Node conformity'):
        errors += Internal.checkPyTree(CTK.t, level=1)
    if (v == 'All conformity' or v == ' > Unique base name'):
        errors += Internal.checkPyTree(CTK.t, level=2)
    if (v == 'All conformity' or v == ' > Unique zone name'):
        errors += Internal.checkPyTree(CTK.t, level=3)
    if (v == 'All conformity' or v == ' > Unique BC name'):
        errors += Internal.checkPyTree(CTK.t, level=4)
    if (v == 'All conformity' or v == ' > BC ranges'):
        errors += Internal.checkPyTree(CTK.t, level=5)
    if (v == 'All conformity' or v == ' > BC match'):
        errors += Internal.checkPyTree(CTK.t, level=6)
    if (v == 'All conformity' or v == ' > Referenced families'):
        errors += Internal.checkPyTree(CTK.t, level=7)
    if (v == 'All conformity' or v == ' > Valid CGNS types'):
        errors += Internal.checkPyTree(CTK.t, level=8)
    if (v == 'All conformity' or v == ' > Valid element nodes'):
        errors += Internal.checkPyTree(CTK.t, level=9)
    if (v == 'All conformity' or v == ' > Valid CGNS flowfield name'):
        errors += Internal.checkPyTree(CTK.t, level=10)
    if (v == 'Multigrid compatibility'):
        MGlevel = CTK.varsFromWidget(VARS[2].get(), type=2)
        minBlk = CTK.varsFromWidget(VARS[0].get(), type=2)
        minBC = CTK.varsFromWidget(VARS[1].get(), type=2)
        if len(MGlevel)>0 and len(minBlk)>0 and len(minBC)>0:
            errors += Internal.checkMultigrid(CTK.t, level=MGlevel[0],
                                              nbMinCoarseB=minBlk[0], nbMinCoarseW=minBC[0])
    if (v == 'Maximum number of nodes'):
        minBlk = CTK.varsFromWidget(VARS[0].get(), type=2)
        if len(minBlk)>0:
            errors = Internal.checkSize(CTK.t, sizeMax=minBlk[0])
        
    if len(errors) == 0: errors = [0, 'No error found.']
    Panels.displayErrors(errors, header='Checking pyTree')
    CTK.TXT.insert('START', 'pyTree checked.\n')

#==============================================================================
def correctPyTree():
    if CTK.t == []: return
    v = VARS[3].get()
    if v == 'Multigrid compatibility':
        CTK.TXT.insert('START', 'Can not correct for multigrid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if v == 'Maximum number of nodes':
        CTK.TXT.insert('START', 'Can not correct for maximum number of nodes.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    if (v == 'All conformity' or v == ' > Node conformity'):
        Internal._correctPyTree(CTK.t, level=1)
    if (v == 'All conformity' or v == ' > Unique base name'):
        Internal._correctPyTree(CTK.t, level=2)
    if (v == 'All conformity' or v == ' > Unique zone name'):
        Internal._correctPyTree(CTK.t, level=3)
    if (v == 'All conformity' or v == ' > Unique BC name'):
        Internal._correctPyTree(CTK.t, level=4)
    if (v == 'All conformity' or v == ' > BC ranges'):
        Internal._correctPyTree(CTK.t, level=5)
    if (v == 'All conformity' or v == ' > BC match'):
        Internal._correctPyTree(CTK.t, level=6)
    if (v == 'All conformity' or v == ' > Referenced families'):
        Internal._correctPyTree(CTK.t, level=7)
    if (v == 'All conformity' or v == ' > Valid CGNS types'):
        Internal._correctPyTree(CTK.t, level=8)
    if (v == 'All conformity' or v == ' > Valid element nodes'):
        Internal._correctPyTree(CTK.t, level=9)
    if (v == 'All conformity' or v == ' > Valid CGNS flowfield name'):
        Internal.correctPyTree(CTK.t, level=10)
    CTK.t = C.fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CTK.TXT.insert('START', 'pyTree corrected.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCheckPyTree', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Check your pyTree.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Button-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=1)
    Frame.columnconfigure(4, weight=1)
    Frame.columnconfigure(5, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkCheckPyTree')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- MinPtsPerCoarseGrid -
    V = TK.StringVar(win); V.set('5'); VARS.append(V)
    if CTK.PREFS.has_key('tkCheckPyTree0'): V.set(CTK.PREFS['tkCheckPyTree0'])
    # -1- MinPtsPerCoarseWin -
    V = TK.StringVar(win); V.set('3'); VARS.append(V)
    if CTK.PREFS.has_key('tkCheckPyTree1'): V.set(CTK.PREFS['tkCheckPyTree1'])
    # -2- Multigrid level
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    if CTK.PREFS.has_key('tkCheckPyTree2'): V.set(CTK.PREFS['tkCheckPyTree2'])
    # -3- global option menu -> things to check
    V = TK.StringVar(win); V.set('All conformity'); VARS.append(V)

    # - MinPtsPerZone on coarse level -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    BB = CTK.infoBulle(parent=B,
                       text='Checked maximum number of points per zone, \nif Maximum number of nodes is selected. \nChecked minimum number of points per direction on coarse grid \nwhen Multigrid compatibility is selected.')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)

    # - MinPtsPerWin on coarse level -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='Checked minimum number of points per direction in BCs on coarse grid, \nwhen Multigrid compatibility is selected.')
    B.grid(row=0, column=2, columnspan=2, sticky=TK.EW)

    # - Checked multigrid level -
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='Checked multigrid level,\nwhen Multigrid compatibility is selected.')
    B.grid(row=0, column=4, columnspan=2, sticky=TK.EW)
    
    # Option menu
    B = TTK.OptionMenu(Frame, VARS[3], 'All conformity', ' > Node conformity', ' > Unique base name', ' > Unique zone name', ' > Unique BC name', 
                       ' > BC ranges', ' > BC match', ' > Referenced families', ' > Valid CGNS types', ' > Valid element nodes',
                       ' > Valid CGNS flowfield name', 'Multigrid compatibility', 'Maximum number of nodes')
    B.grid(row=1, column=0, columnspan=8, sticky=TK.EW)

    # - Buttons -
    B = TTK.Button(Frame, text="Check", command=runCheckPyTree)
    B.grid(row=2, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Check pyTree.')
    B = TTK.Button(Frame, text="Correct", command=correctPyTree)
    B.grid(row=2, column=4, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Attemp to correct pyTree.')
    
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
def resetApp():
    VARS[0].set('5')
    VARS[1].set('3')
    VARS[2].set('1')
    CTK.PREFS['tkCheckPyTree0'] = VARS[0].get()
    CTK.PREFS['tkCheckPyTree1'] = VARS[1].get()
    CTK.PREFS['tkCheckPyTree2'] = VARS[2].get()
    savePrefFile()

#==============================================================================
def saveApp():
    CTK.PREFS['tkCheckPyTree0'] = VARS[0].get()
    CTK.PREFS['tkCheckPyTree1'] = VARS[1].get()
    CTK.PREFS['tkCheckPyTree2'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkCheckPyTree '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
