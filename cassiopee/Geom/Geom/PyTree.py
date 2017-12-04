"""Geometry definition module.
"""
# 
# Python Interface to define geometries in pyTrees
#
import Geom
__version__ = Geom.__version__

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
except:
    raise ImportError("Geom.PyTree: requires Converter.PyTree module.")

def point(P):
    """Create a point.
    Usage: point((x,y,z))"""
    a = Geom.point(P)
    return C.convertArrays2ZoneNode('point', [a])

def naca(epaisseur, N=101):
    """Create a naca profile of N points.
    Usage: naca(epaisseur, N)"""
    a = Geom.naca(epaisseur, N)
    return C.convertArrays2ZoneNode('naca', [a])    
    
def line(P1, P2, N=100):
    """Create a line of N points. Usage: line( (x1,y1,z1), (x2,y2,z2), N )"""
    a = Geom.line(P1, P2, N)
    return C.convertArrays2ZoneNode('line', [a])


def polyline(Pts):
    """Create a polyline of N points.
    Usage: polyline( [(x1,y1,z1),....,(xn,yn,zn)])"""
    a = Geom.polyline( Pts )
    return C.convertArrays2ZoneNode('polyline', [a])

def circle(Center, R, tetas=0., tetae=360., N=100):
    """Create a portion of circle of N points and of center C,
    radius R, between angle tetas and tetae.
    Usage: circle((xc,yc,zc), R, tetas, tetae, N)"""
    a = Geom.circle(Center, R, tetas, tetae, N)
    return C.convertArrays2ZoneNode('circle', [a])

def bezier(t, N=100, M=100, density=-1):
    """Create a a Bezier curve controlled by an array of control points.
    Usage: bezier(tc, N)"""
    return C.TZGC(t, 'nodes', Geom.bezier, N, M, density)

def spline(t, order=3, N=100, M=100, density=-1):
    """Create a spline of N points.
    Usage: spline(ctrlsPts, order, N)"""
    return C.TZGC(t, 'nodes', Geom.spline, order, N, M, density)

def nurbs(t, weight='weight', order=3, N=100, M=100, density=-1):
    """Create a nurbs of N points.
    Usage: nurbs(ctrlPts, order, N)"""
    w = C.getField(weight,t)[0]
    Pts = C.getFields(Internal.__GridCoordinates__, t)[0]
    Pts = Converter.addVars([Pts,w])
    surf = Geom.nurbs(Pts, weight, order, N, M, density)
    return C.convertArrays2ZoneNode('nurbs', [surf])
    
def curve(f, N=100):
    """Create a curve from a user defined parametric function.
    Usage: curve(f, N)"""
    a = Geom.curve(f, N)
    return C.convertArrays2ZoneNode('curve', [a])

def cone(center, Rb, Rv, H, N=100):
    """Create cone of NxNh points and of center C, basis radius Rb,
    vertex radius Rv and height H.
    Usage: cone((xc,yc,zc), Rb, Rv, H, N)"""
    a = Geom.cone(center, Rb, Rv, H, N)
    return C.convertArrays2ZoneNode('cone', [a])

def sphere(center, R, N=100):
    """Create a sphere of Nx2N points and of center C and radius R.
    Usage: sphere((xc,yc,zc), R, N)"""
    a = Geom.sphere(center, R, N)
    return C.convertArrays2ZoneNode('sphere', [a])

def sphere6(center, R, N=100):
    """Create a shpere of NxN points and of center C and radius R.
    Usage: sphere((xc,yc,zc), R, N)"""
    A = Geom.sphere6(center, R, N)
    return [C.convertArrays2ZoneNode('sphere-part1', [A[0]]),
            C.convertArrays2ZoneNode('sphere-part2', [A[1]]),
            C.convertArrays2ZoneNode('sphere-part3', [A[2]]),
            C.convertArrays2ZoneNode('sphere-part4', [A[3]]),
            C.convertArrays2ZoneNode('sphere-part5', [A[4]]),
            C.convertArrays2ZoneNode('sphere-part6', [A[5]])]

def sphereYinYang(center, R, N=100):
    """Create a sphere of center C and radius R made of two overlapping zones.
    Usage: sphereYinYang((xc,yc,zc), R, N)"""
    A = Geom.sphereYinYang(center, R, N)
    return [C.convertArrays2ZoneNode('sphere-part1', [A[0]]),
            C.convertArrays2ZoneNode('sphere-part2', [A[1]])]

def torus(center, R, r, alphas=0., alphae=360.,
          betas=0., betae=360., NR=100, Nr=100):
    """Create NRxNr points lying on a torus of center C and radii R (main)
    and r (tube) between the angles alphas and alphae (XY-plane) and
    between betas and betae (RZ-plane).
    Usage: torus((xc,yc,zc), R, r, NR, Nr, alphas, alphae, betas, betae)"""
    a = Geom.torus(center, R, r, alphas, alphae, betas, betae, NR, Nr)
    return C.convertArrays2ZoneNode('torus', [a])
            
def triangle(P1, P2, P3):
    """Create a single triangle with points P1, P2, P3.
    Usage: triangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3))"""
    a = Geom.triangle(P1, P2, P3)
    return C.convertArrays2ZoneNode('triangle', [a])

def quadrangle(P1, P2, P3, P4):
    """Create a single quadrangle with points P1, P2, P3, P4.
    Usage: quadrangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4))"""
    a = Geom.quadrangle(P1, P2, P3, P4)
    return C.convertArrays2ZoneNode('quadrangle', [a])

def surface(f, N=100):
    """Create a surface from a user defined parametric function.
    Usage: surface(f, N)"""
    a = Geom.surface(f, N)
    return C.convertArrays2ZoneNode('surface', [a])

def getLength(t):
    """Return the length of 1D array(s) defining a mesh.
    Usage: getLength(t)"""
    coords = C.getFields(Internal.__GridCoordinates__, t)
    return Geom.getLength(coords)

def getDistantIndex(t, ind, l):
    """Return the index of 1D array defining a mesh located at a
    distance l of ind.
    Usage: getDistantIndex(t, ind, l)"""
    a = C.getFields(Internal.__GridCoordinates__, t)[0]
    return Geom.getDistantIndex(a, ind, l)

def getNearestPointIndex(t, pointList):
    """Return the nearest index of points in array.
    Usage: getNearestPointIndex(t, pointList)"""
    a = C.getFields(Internal.__GridCoordinates__, t)[0]
    return Geom.getNearestPointIndex( a, pointList )

def getCurvatureHeight(t):
    """Return the curvature height for each point.
    Usage: getCurvatureHeight(t)"""
    return C.TZGC(t, 'nodes', Geom.getCurvatureHeight)

def _getCurvatureHeight(t):
    return C._TZGC(t, 'nodes', Geom.getCurvatureHeight)

def getCurvatureRadius(t):
    """Return the curvature radius for each point.
    Usage: getCurvatureRadius(t)"""
    return C.TZGC(t, 'nodes', Geom.getCurvatureRadius)
    
def _getCurvatureRadius(t):
    return C._TZGC(t, 'nodes', Geom.getCurvatureRadius)

def getCurvatureAngle(t):
    """Return the curvature angle for each point...
    Usage: getCurvatureAngle(t)"""
    return C.TZGC(t, 'nodes', Geom.getCurvatureAngle)

def _getCurvatureAngle(t):
    return C._TZGC(t, 'nodes', Geom.getCurvatureAngle)

def getSharpestAngle(t):
    """Return the sharpest angle for each point of a surface based on the sharpest angle
    between adjacent element to which the point belongs to.
    Usage: getSharpestAngle(a)"""
    return C.TZGC(t, 'nodes', Geom.getSharpestAngle)

def _getSharpestAngle(t):
    return C._TZGC(t, 'nodes', Geom.getSharpestAngle)

def getCurvilinearAbscissa(t):
    """Return the curvilinear abscissa for each point...
    Usage: getCurvilinearAbscissa(t)"""
    return C.TZGC(t, 'nodes', Geom.getCurvilinearAbscissa)
    
def getDistribution(t):
    """Return the curvilinear abscissa for each point as coordinates
    Usage: getDistribution(t)"""
    return C.TZGC(t, 'nodes', Geom.getDistribution)    

def _getCurvilinearAbscissa(t):
    return C._TZGC(t, 'nodes', Geom.getCurvilinearAbscissa)

def getTangent(t):
    """
    Makes the tangent of a 1D curve. The input argument shall be a structured
    1D curve. Each node of the output represents the unitary tangent vector, 
    pointing towards the tangent direction of the input 1D curve.
    Usage: b = getTangent(t)"""
    tp = Internal.copyRef(t)
    C._deleteFlowSolutions__(tp)
    C._TZGC(tp, 'nodes', Geom.getTangent)
    return tp

def addSeparationLine(t, line0):
    """Add a separation line defined in line0 to a mesh defined in t.
    Usage: addSeparationLine(t, line0)"""
    al = C.getFields(Internal.__GridCoordinates__, line0)[0]
    at = C.getFields(Internal.__GridCoordinates__, t)[0]
    arrays = Geom.addSeparationLine(at, al)
    zones = []
    for i in arrays:
        zone = C.convertArrays2ZoneNode(t[0], [i])        
        zones.append(zone)
    return zones

def lineGenerate(t, line):
    """Generate a surface mesh by using 1D array (defining a mesh)
    and following the curve defined in line.
    Usage: lineGenerate(t, line)"""
    al = C.getFields(Internal.__GridCoordinates__, line)
    if len(al) == 1: al = al[0]
    al2 = Converter.node2Center(al)
    # Attention les coord. des centres ne sont pas justes! mais
    # elles ne sont pas utilisees dans la fonction
    return C.TZAGC(t, 'both', 'both', Geom.lineGenerate,
                   Geom.lineGenerate, al, al2)

def axisym(t, center, axis, angle=360., Ntheta=360, rmod=None):
    """Create an axisymmetric mesh given an azimuthal surface mesh.
    Usage: axisym(t, (xo,yo,zo), (nx,ny,nz), teta, Nteta, rmod)"""
    # Attention en centres, les coord. des centres ne sont pas justes! mais
    # elles ne sont pas utilisees dans la fonction
    if rmod is not None: rmod = C.getFields(Internal.__GridCoordinates__, rmod)[0]
    return C.TZAGC(t, 'both', 'both', Geom.axisym, Geom.axisym,
                   center, axis, angle, Ntheta, rmod,
                   center, axis, angle, Ntheta-1, rmod)

def _axisym(t, center, axis, angle=360., Ntheta=360, rmod=None):
    if rmod is not None: rmod = C.getFields(Internal.__GridCoordinates__, rmod)[0]
    return C._TZAGC(t, 'both', 'both', Geom.axisym, Geom.axisym,
                    center, axis, angle, Ntheta, rmod,
                    center, axis, angle, Ntheta-1, rmod)

def volumeFromCrossSections(t):
    """Generate a 3D volume from cross sections contours in (x,y) planes.
    Usage: volumeFromCrossSections(t)"""
    tp = Internal.copyRef(t)
    nodes = Internal.getZones(tp)
    coords = []
    for z in nodes:
        coords.append(C.getFields(Internal.__GridCoordinates__, z)[0])
    coordp = Geom.volumeFromCrossSections(coords)
    zone = C.convertArrays2ZoneNode('Volume', [coordp])
    return zone

def text1D(string, font='text1', smooth=0, offset=0.5):
    """Create a 1D text. offset is the space between letters,
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text1D(string, font, smooth, offset)"""
    a = Geom.text1D(string, font, smooth, offset)
    zones = []
    for i in a:
        zone = C.convertArrays2ZoneNode(string, [i])
        zones.append(zone)
    return zones

def text2D(string, font='text1', smooth=0, offset=0.5):
    """Create a 2D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text2D(string, font, smooth, offset)"""
    a = Geom.text2D(string, font, smooth, offset)
    return C.convertArrays2ZoneNode(string, [a])

def text3D(string, font='text1', smooth=0, offset=0.5):
    """Create a 3D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text3D(string, font, smooth, offset)"""
    a = Geom.text3D(string, font, smooth, offset)
    return C.convertArrays2ZoneNode(string, [a])

def connect1D(curves, sharpness=0, N=10, lengthFactor=1.):
    """Connect curves with sharp or smooth junctions.
    Usage: a = connect1D(A, sharpness=0)"""
    a = C.getFields(Internal.__GridCoordinates__, curves)
    z = Geom.connect1D(a, sharpness, N, lengthFactor)
    return C.convertArrays2ZoneNode('connected', [z])
