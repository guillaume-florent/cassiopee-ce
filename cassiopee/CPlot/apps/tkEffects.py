# - render effects (shadow, dof,...) -
import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setShadow(event=None):
    if VARS[1].get() == '1': CPlot.setState(shadow=1)
    else: CPlot.setState(shadow=0)

#==============================================================================
def setDOF(event=None):
    if VARS[2].get() == '1': CPlot.setState(dof=1)
    else: CPlot.setState(dof=0)

#==============================================================================
def setDofPower(event=None):
    #off = WIDGETS['dofPower'].get()*8./100.+2.
    off = WIDGETS['dofPower'].get()*2./100.+2.
    VARS[5].set('Depth of field power [%.2f].'%off)
    CPlot.setState(dofPower=off)
    
#==============================================================================
def setLightOffsetX(event=None):
    off = WIDGETS['lightOffsetX'].get() / 50. - 1.
    VARS[3].set('Light offset in x [%.2f %%].'%off)
    CPlot.setState(lightOffset=(off, -999))

#==============================================================================
def setLightOffsetY(event=None):
    off = WIDGETS['lightOffsetY'].get() / 50.
    VARS[4].set('Light offset in y [%.2f %%].'%off)
    CPlot.setState(lightOffset=(-999, off))
    
#==============================================================================
def setViewAngle(event=None):
    angle = VARS[0].get()
    angle = float(angle)
    CPlot.setState(viewAngle=angle)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkEffects', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set anaglyph mode.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<Button-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkEffects')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- cam view angle -
    V = TK.StringVar(win); V.set('50.'); VARS.append(V)
    if CTK.PREFS.has_key('tkEffectsAngle'): V.set(CTK.PREFS['tkEffectsAngle'])
    # -1- Shadow
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if CTK.PREFS.has_key('tkEffectsShadow'):
        V.set(CTK.PREFS['tkEffectsShadow'])
    # -2- DOF
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if CTK.PREFS.has_key('tkEffectsDOF'): V.set(CTK.PREFS['tkEffectsDOF'])
    # -3- Light offset X info bulle
    V = TK.StringVar(win); V.set('Light offset in x.'); VARS.append(V)
    # -4- Light offset Y info bulle
    V = TK.StringVar(win); V.set('Light offset in y.'); VARS.append(V)
    # -5- Depth of field info bulle
    V = TK.StringVar(win); V.set('Depth of field power.'); VARS.append(V)

    # - Camera angle -
    B = TTK.Button(Frame, text="Set Cam angle", command=setViewAngle)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the camera view angle.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', setViewAngle)
    BB = CTK.infoBulle(parent=B, text='Cam view angle (deg).')

    # - Light offset X -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setLightOffsetX, value=55)
    WIDGETS['lightOffsetX'] = B
    B.grid(row=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[3])
    
    # - Light offset Y -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setLightOffsetY, value=25)
    WIDGETS['lightOffsetY'] = B
    B.grid(row=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[4])
    
    # - DOF power -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setDofPower, value=50)
    WIDGETS['dofPower'] = B
    B.grid(row=3, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[5])
    
    # - Shadow + DOF -
    B = TTK.Checkbutton(Frame, text='Shadow', variable=VARS[1],
                        command=setShadow)
    BB = CTK.infoBulle(parent=B, text='Toggle shadows.')
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Checkbutton(Frame, text='DOF', variable=VARS[2], command=setDOF)
    BB = CTK.infoBulle(parent=B, text='Toggle depth of field blur.')
    B.grid(row=4, column=1, sticky=TK.EW)
    
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
    CTK.PREFS['tkEffectsAngle'] = VARS[0].get()
    CTK.PREFS['tkEffectsShadow'] = VARS[1].get()
    CTK.PREFS['tkEffectsDOF'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('50.')
    VARS[1].set('0')
    VARS[2].set('0')
    CTK.PREFS['tkEffectsAngle'] = VARS[0].get()
    CTK.PREFS['tkEffectsShadow'] = VARS[1].get()
    CTK.PREFS['tkEffectsDOF'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkStereo '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
