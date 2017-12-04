"""Plotter functions.
"""
__version__ = '2.5'
__author__ = "Christophe Benoit, Stephanie Peron, Pascal Raud, Matthieu Soismier, Bertrand Michel"
#
# Plotter for arrays
#
import cplot
import time
__timeStep__ = 0.02
__slot__ = None

#==============================================================================
def configure(use_render):
    """
    Configure CPlot pour soit un rendering direct (cplot.useDirect),
    utiliser les Display Lists (cplot.useDL)
    soit les VBO (cplt.useVBO)
    """
    # Plus tard, cplt.useImg egalement ?
    cplot.configure(use_render)

#==============================================================================
# -- display --
def display(arrays,
            dim=-1,
            mode=-1,
            scalarField=-1,
            vectorField1=-1, vectorField2=-1, vectorField3=-1,
            displayBB=-1,
            displayInfo=-1,
            displayIsoLegend=-1,
            meshStyle=-1,
            solidStyle=-1,
            scalarStyle=-1,
            vectorStyle=-1,
            colormap=-1,
            niso=-1,
            isoEdges=-1,
            isoScales=[],
            win=(-1,-1),
            posCam=(-999,-999,-999),
            posEye=(-999,-999,-999),
            dirCam=(-999,-999,-999),
            viewAngle=-1.,
            bgColor=-1,
            shadow=-1,
            dof=-1,
            stereo=-1,
            stereoDist=-1.,
            export="None",
            exportResolution="None",
            zoneNames=[],
            renderTags=[],
            offscreen=0):
    """Display arrays.
    Usage: display(arrays)"""
    if arrays != [] and not isinstance(arrays[0], list): arrays = [arrays]
    if __slot__ is None:
        import os.path
        shaderPath = os.path.dirname(__file__)+'/'
        cplot.setShaderPath(shaderPath)
        displayNew__(arrays, dim, mode, scalarField, vectorField1,
                     vectorField2, vectorField3, displayBB, displayInfo,
                     displayIsoLegend, meshStyle, solidStyle,
                     scalarStyle, vectorStyle, colormap,
                     niso, isoEdges, isoScales, win,
                     posCam, posEye, dirCam, viewAngle, bgColor,
                     shadow, dof, stereo, stereoDist,
                     export, exportResolution, zoneNames, renderTags,
                     offscreen)
    else:
        displayAgain__(arrays, dim, mode, scalarField, vectorField1,
                       vectorField2, vectorField3, displayBB, displayInfo,
                       displayIsoLegend, meshStyle, solidStyle,
                       scalarStyle, vectorStyle,
                       colormap, niso, isoEdges, isoScales, win,
                       posCam, posEye, dirCam, viewAngle, bgColor,
                       shadow, dof, stereo, stereoDist,
                       export, exportResolution,
                       zoneNames, renderTags, offscreen)

#==============================================================================
def render():
    """Force render.
    Usage: render()"""
    cplot.render()

#==============================================================================
def delete(list):
    """Delete zones from plotter.
    Usage: delete([i1,i2,...])"""
    if __slot__ is None: return
    cplot.delete(list)

#==============================================================================
def add(arrays, no, array, zoneName=None, renderTag=None):
    """Add one zone to plotter.
    Usage: add(arrays, no, array, zoneName, renderTag)"""
    nzs = 0; nzu = 0
    for i in xrange(no):
        if len(arrays[i]) == 5: nzs += 1
        else: nzu += 1
    arrays.insert(no, array)
    if __slot__ is None:
        if (no == -1): no = len(arrays)
        arrays.insert(no, array)
        if len(arrays) == 1:
            display(arrays, zoneNames=[zoneName], renderTags=[renderTag])
        else: display(arrays)
    else:
        cplot.add(array, (nzs, nzu), zoneName, renderTag)

#==============================================================================
def replace(arrays, no, array, zoneName=None, renderTag=None):
    """Replace arrays[no] by array.
    Usage: replace(arrays, no, array, zoneName, renderTag)"""
    zone = arrays[no]
    if len(zone) == 5: oldType = 1
    else: oldType = 2
    nzs = 0; nzu = 0
    for i in xrange(no):
        if (len(arrays[i]) == 5): nzs += 1
        else: nzu += 1
    arrays[no] = array
    if __slot__ is None:
        if len(arrays) == 1:
            display(arrays, zoneNames=[zoneName], renderTags=[renderTag])
        else: display(arrays)
    else:
        cplot.replace(array, (nzs, nzu, oldType), zoneName, renderTag)

#==============================================================================
def display1D(arrays, slot=0, gridPos=(0,0), gridSize=(-1,-1),
              bgBlend=1., var1='x', var2='y', r1=None, r2=None):
    """Display 1D plots.
    Usage: display1D(arrays, slot, ....)"""
    import Converter
    import numpy
    if len(arrays) > 0:
        if isinstance(arrays[0], numpy.ndarray): # numpy [x,y]
            if len(arrays) < 2:
                raise ValueError, 'display1D: requires at least two numpys [x,y]'
            x = arrays[0]; y = arrays[1]
            n = x.size
            array = Converter.array(var1+','+var2,n,1,1)
            array[1][0,:] = x[:]
            array[1][1,:] = y[:]
            arrays = [array]
        elif not isinstance(arrays[0], list): arrays = [arrays]
    if __slot__ is None: display([])
    try:
        arrays = Converter.convertArray2Hexa(arrays)
        minv1 = Converter.getMinValue(arrays, var1)
        maxv1 = Converter.getMaxValue(arrays, var1)
        minv2 = Converter.getMinValue(arrays, var2)
        maxv2 = Converter.getMaxValue(arrays, var2)
        if r1 is None: r1 = (minv1, maxv1)
        if r2 is None: r2 = (minv2, maxv2)
    except: pass
    if r1 is None: r1 = (0.,1.)
    if r2 is None: r2 = (0.,1.)
    cplot.display1D(arrays, slot, gridPos, gridSize,
                    bgBlend, var1, var2, r1, r2)
    time.sleep(__timeStep__)

#==============================================================================
def pressKey():
    """Wait for user to press a key.
    Usage: pressKey()"""
    cplot.pressKey()

#==============================================================================
# -- get/set functions --
def getState(mode):
    """Return a state in plotter.
    Usage: n = getState(mode)"""
    return cplot.getState(mode)

def getSelectedZone():
    """Return the selected zone in plotter.
    Usage: n = getSelectedZone()"""
    return cplot.getSelectedZone()

def getSelectedZones():
    """Return the selected zones in plotter.
    Usage: list = getSelectedZones()"""
    return cplot.getSelectedZones()

def getSelectedStatus(zone):
    """Return the selected status of a zone in plotter.
    Usage: status = getSelectedStatus(zone)"""
    return cplot.getSelectedStatus(zone)

def getActiveZones():
    """Return the active (displayed) zones in plotter.
    Usage: list = getActiveZones()"""
    return cplot.getActiveZones()

def getActiveStatus(zone):
    """Return the active status of a zone in plotter.
    Usage: status = getActiveStatus(zone)"""
    return cplot.getActiveStatus(zone)

def getActivePoint():
    """Return the active (clicked) point in plotter.
    Usage: n = getActivePoint()"""
    return cplot.getActivePoint()

def getActivePointIndex():
    """Return the active (clicked) point index.
    Usage: n = getActivePointIndex()"""
    return cplot.getActivePointIndex()

def getActivePointF():
    """Return the active (clicked) point field values.
    Usage: f = getActivePointF()"""
    return cplot.getActivePointF()

def getMouseState():
    """Return mouse state (mouse position and button state)."""
    return cplot.getMouseState()

def getKeyboard():
    """Return the pressed keys.
    Usage: n = getKeyboard()"""
    return cplot.getKeyboard()

def resetKeyboard():
    """Reset the keyboard string.
    Usage: resetKeyboard()"""
    return cplot.resetKeyboard()

def setState(dim=-1,
             mode=-1,
             scalarField=-1,
             vectorField1=-1, vectorField2=-1, vectorField3=-1,
             displayBB=-1,
             displayInfo=-1,
             displayIsoLegend=-1,
             meshStyle=-1,
             solidStyle=-1,
             scalarStyle=-1,
             vectorStyle=-1,
             colormap=-1,
             niso=-1,
             isoEdges=-1,
             isoScales=[],
             win=(-1,-1),
             posCam=(-999,-999,-999),
             posEye=(-999,-999,-999),
             dirCam=(-999,-999,-999),
             viewAngle=-1,
             lightOffset=(-999,-999),
             bgColor=-1,
             shadow=-1,
             dof=-1, dofPower=-1,
             ghostifyDeactivatedZones=-1,
             edgifyActivatedZones=-1,
             edgifyDeactivatedZones=-1,
             export="None",
             exportResolution="None",
             continuousExport=-1,
             envmap="None", message="None",
             stereo=-1, stereoDist=-1.,
             cursor=-1,
             gridSize=(-1, -1),
             timer=-1,
             selectionStyle=-1,
             activateShortCuts=-1,
             billBoards=None):
    """Set CPlot state.
    Usage: setState(posCam=(12,0,0))"""
    cplot.setState(dim, mode, scalarField, vectorField1, vectorField2,
                   vectorField3, displayBB, displayInfo, displayIsoLegend,
                   meshStyle, solidStyle, scalarStyle,
                   vectorStyle, colormap,
                   niso, isoEdges, isoScales, win,
                   posCam, posEye, dirCam, viewAngle, lightOffset,
                   bgColor, shadow, dof, dofPower,
                   ghostifyDeactivatedZones, edgifyActivatedZones,
                   edgifyDeactivatedZones,
                   export, exportResolution, continuousExport,
                   envmap, message,
                   stereo, stereoDist, cursor, gridSize, timer, selectionStyle,
                   activateShortCuts, billBoards)

def setMode(mode):
    """Set CPlot display mode.
    Usage: setMode(0)"""
    cplot.setMode(mode)

def changeVariable():
    """Change displayed variable.
    Usage: changeVariable()"""
    cplot.changeVariable()

def changeStyle():
    """Change CPlot display style.
    Usage: changeStyle()"""
    cplot.changeStyle()

def changeInfoDisplay():
    """Change CPlot info display style.
    Usage: changeInfoDisplay()"""
    cplot.changeInfoDisplay()

def changeBlanking():
    """Change the blanking procedure.
    Usage: changeBlanking()"""
    cplot.changeBlanking()

def setDim(dim):
    """Set CPlot display dim 3, 2 or 1.
    Usage: setDim(2)"""
    cplot.setDim(dim)

def setActivePoint(x,y,z):
    """Set the active (clicked) point in plotter.
    Usage: setActivePoint(x,y,z)"""
    return cplot.setActivePoint(x,y,z)

def setSelectedZones(list):
    """Set selected zones.
    Usage: setSelectedZones([(0,1),(1,1),...])"""
    cplot.setSelectedZones(list)

def unselectAllZones():
    """Unselect all zones.
    Usage: unselectAllZones()"""
    cplot.unselectAllZones()

def setActiveZones(list):
    """Set active (displayed) zones.
    Usage: setActiveZones([(0,1)])"""
    cplot.setActiveZones(list)

def setZoneNames(list):
    """Set zone names.
    Usage: setZoneNames([(0,'myZone')])"""
    cplot.setZoneNames(list)

#==============================================================================
def lookFor():
    """Look for selected zones.
    Usage: lookFor()"""
    cplot.lookFor()

def fitView():
    """Fit the view to objects.
    Usage: fitView()"""
    cplot.fitView()

def finalizeExport(t=0):
    """Finalize export for continuous export."""
    cplot.finalizeExport(t)

def hide():
    """Hide window."""
    cplot.hide()

def show():
    """Show window if it has been hidden with flush."""
    cplot.show()

#==============================================================================
# camera
#==============================================================================
def moveCamera(checkPoints, moveEye=False, N=100, speed=0.001):
    import Geom
    N = max(N, 3)
    p = Geom.polyline(checkPoints)
    if len(checkPoints) == 2: d = Geom.spline(p, 2, N)
    else: d = Geom.spline(p, 3, N)
    i = 0
    while i < N-1:
        time.sleep(__timeStep__*speed*0.06)
        if i > N-11: inc = N-i-1
        else: inc = 10
        posCam = (d[1][0,i],d[1][1,i],d[1][2,i])
        if moveEye: 
            posEye = (d[1][0,i+inc],d[1][1,i+inc],d[1][2,i+inc])
            setState(posCam=posCam, posEye=posEye)
        else: setState(posCam=posCam)
        i += 1

def travelRight(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2.*3.14*R*xr*0.5
    d2 = Vector.cross(d1, dirCam)
    d2 = Vector.normalize(d2)
    d2 = Vector.mul(L,d2)
    P2 = Vector.sub(posCam, d2)
    d3 = Vector.sub(posEye, P2)
    d4 = Vector.cross(d3, dirCam)
    d4 = Vector.normalize(d4)
    d4 = Vector.mul(L,d4)
    P3 = Vector.sub(P2, d4)
    checkPoints = [posCam,tuple(P2),tuple(P3)]
    moveCamera(checkPoints, N=N)

def travelLeft(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*0.5
    d2 = Vector.cross(d1, dirCam)
    d2 = Vector.normalize(d2)
    d2 = Vector.mul(L,d2)
    P2 = Vector.add(posCam, d2)
    d3 = Vector.sub(posEye, P2)
    d4 = Vector.cross(d3, dirCam)
    d4 = Vector.normalize(d4)
    d4 = Vector.mul(L,d4)
    P3 = Vector.add(P2, d4)
    checkPoints = [posCam,tuple(P2),tuple(P3)]
    moveCamera(checkPoints, N=N)

def travelUp(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*1.
    d2 = Vector.normalize(dirCam)
    d2 = Vector.mul(L,d2)
    P2 = Vector.add(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    moveCamera(checkPoints, N=N)

def travelDown(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*1.
    d2 = Vector.normalize(dirCam)
    d2 = Vector.mul(L,d2)
    P2 = Vector.sub(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    moveCamera(checkPoints, N=N)

def travelIn(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = R*xr
    d2 = Vector.mul(L,d1)
    P2 = Vector.add(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    moveCamera(checkPoints, N=N)

def travelOut(xr=0.1, N=100):
    import KCore.Vector as Vector
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = R*xr
    d2 = Vector.mul(L,d1)
    P2 = Vector.sub(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    moveCamera(checkPoints, N=N)

#==============================================================================

# -- Internal functions --
def setFileName__(name):
    cplot.setFileName(name)

def displayNew__(arrays, dim, mode, scalarField, vectorField1, vectorField2,
                 vectorField3, displayBB, displayInfo, displayIsoLegend,
                 meshStyle, solidStyle, scalarStyle, vectorStyle,
                 colormap, niso, isoEdges, isoScales, win,
                 posCam, posEye, dirCam, viewAngle, bgColor,
                 shadow, dof, stereo, stereoDist,
                 export, exportResolution, zoneNames, renderTags, offscreen):
    global __slot__
    import threading
    a = threading.Thread(None, cplot.displayNew, None,
                         (arrays, dim, mode, scalarField, vectorField1,
                          vectorField2, vectorField3, displayBB, displayInfo,
                          displayIsoLegend,
                          meshStyle, solidStyle, scalarStyle,
                          vectorStyle, colormap, niso, isoEdges, isoScales,
                          win, posCam, posEye, dirCam, viewAngle, bgColor,
                          shadow, dof, stereo, stereoDist,
                          export, exportResolution,
                          zoneNames, renderTags, offscreen), {})
    a.start()
    __slot__ = a

def displayAgain__(arrays, dim, mode, scalarField, vectorField1, vectorField2,
                   vectorField3, displayBB, displayInfo, displayIsoLegend,
                   meshStyle, solidStyle, scalarStyle, vectorStyle,
                   colormap, niso, isoEdges, isoScales,
                   win, posCam, posEye, dirCam, viewAngle, bgColor,
                   shadow, dof, stereo, stereoDist,
                   export, exportResolution, zoneNames, renderTags, offscreen):
    cplot.displayAgain(arrays, dim, mode, scalarField, vectorField1,
                       vectorField2, vectorField3, displayBB, displayInfo,
                       displayIsoLegend,
                       meshStyle, solidStyle, scalarStyle, vectorStyle,
                       colormap, niso, isoEdges, isoScales,
                       win, posCam, posEye, dirCam, viewAngle, bgColor,
                       shadow, dof, stereo, stereoDist,
                       export, exportResolution, zoneNames, renderTags,
                       offscreen)
    time.sleep(__timeStep__)
